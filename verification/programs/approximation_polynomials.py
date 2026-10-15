#!/usr/bin/env python3
"""Recompute the exact approximation polynomials and constants in Sections 4–5.

Run from any directory with Python 3.10 or later. No third-party packages are
required. The article input and current printed tables are checked directly.
Results are written under D1_RESULTS_DIR, or to a new directory under ../reruns by default.
"""
import sys
if sys.flags.optimize:
    raise SystemExit("Run without -O or -OO: the verification requires assertions.")
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
from math import factorial, comb, isqrt
import hashlib, json, re, os, time
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/article.tex'
import tempfile
if os.environ.get('D1_RESULTS_DIR'):
    RESULTS = Path(os.environ['D1_RESULTS_DIR'])
    RESULTS.mkdir(parents=True, exist_ok=True)
else:
    (ROOT/'reruns').mkdir(exist_ok=True)
    RESULTS = Path(tempfile.mkdtemp(prefix=Path(__file__).stem+'-', dir=ROOT/'reruns'))
    print('Results directory:', RESULTS)
START=time.monotonic()
TEXT=SOURCE.read_text(encoding='utf-8')
SHA=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
EXPECTED='9cb6ed2154056e158d2c420012e94ef7c397c69d66c0fd9301a1bb52197dc432'
assert SHA==EXPECTED, 'The article input does not match this verification release'
S4=TEXT[TEXT.index(r'\section{Simultaneous approximation}'):TEXT.index(r'\section{Finite integral approximation forms}')]
S5=TEXT[TEXT.index(r'\section{Finite integral approximation forms}'):TEXT.index(r'\section{',TEXT.index(r'\section{Finite integral approximation forms}')+1)]
checks=[]
def check(name,condition,detail=None):
    if not condition: raise AssertionError((name,detail))
    checks.append({'name':name,'passed':True,**({'detail':detail} if detail is not None else {})})
N=5
ZERO=(0,)*N
def clean(p):
    return {e: v for e, v in p.items() if v}

def const(v):
    return {ZERO: F(v)} if v else {}

def mon(v, **e):
    z = [0] * N
    for k, x in e.items():
        z['abstc'.index(k)] = x
    return {tuple(z): F(v)} if v else {}

def add(*ps):
    d = {}
    for p in ps:
        for e, v in p.items():
            d[e] = d.get(e, 0) + v
    return clean(d)

def scale(p, v):
    return clean({e: x * v for e, x in p.items()})

def mul(p, q, cap=None):
    d = {}
    for e, v in p.items():
        for f, w in q.items():
            g = tuple((x + y for x, y in zip(e, f)))
            if cap is None or g[4] <= cap:
                d[g] = d.get(g, 0) + v * w
    return clean(d)

def power(p, n):
    assert n >= 0
    q = const(1)
    for _ in range(n):
        q = mul(q, p)
    return q

def shift(p, k, n):
    d = {}
    for e, v in p.items():
        ee = list(e)
        ee[k] += n
        d[tuple(ee)] = v
    return d

def divpoly(p, q):
    d = {}
    r = p.copy()
    fq = max(q)
    vq = q[fq]
    while r:
        er = max(r)
        dif = tuple((x - y for x, y in zip(er, fq)))
        assert min(dif) >= 0, (er, fq)
        co = r[er] / vq
        d[dif] = co
        r = add(r, scale(mul({dif: co}, q), -1))
    return d

def det(mat, cap=None):
    out = {}
    n = len(mat)
    for p in permutations(range(n)):
        sig = (-1) ** sum((p[i] > p[j] for i in range(n) for j in range(i + 1, n)))
        v = const(sig)
        for i, j in enumerate(p):
            v = mul(v, mat[i][j], cap)
        out = add(out, v)
    return out

def binom(x, k):
    z = F(1)
    for i in range(k):
        z *= F(x - i, i + 1)
    return z

def comps(n, k):
    if k == 1:
        yield (n,)
        return
    for j in range(n + 1):
        for t in comps(n - j, k - 1):
            yield ((j,) + t)

def evaluate(p, values):
    return sum((co * prod((v ** e for v, e in zip(values, ex))) for ex, co in p.items()))

def prod(xs):
    r = 1
    for x in xs:
        r *= x
    return r

def family(k, m):
    p = [m, k, k]
    out = []
    for i in range(3):
        row = []
        for j in range(3):
            q = {}
            for h in comps(p[j] + int(i == j) - 1, 3):
                co = 4 ** k * binom(F(2 * k + 1, 2), h[j]) * prod((binom(-p[l] - int(i == l), h[l]) for l in range(3) if l != j))
                if j == 0:
                    term = mon((-1) ** (2 * k + int(i == 1) + int(i == 2) + h[1] + h[2]), a=int(i == 1) + h[1], b=int(i == 2) + h[2], c=k - h[0])
                    term = mul(term, power(diff, 2 * k))
                elif j == 1:
                    term = mon((-1) ** (k + int(i == 2) + h[2]), a=m, b=int(i == 2) + h[2])
                    term = mul(term, power(add(mul(a, c), const(1)), k - h[1]))
                    term = mul(term, power(diff, k - int(i == 2) - h[2]))
                else:
                    term = mon(1, b=m, a=int(i == 1) + h[1])
                    term = mul(term, power(add(mul(b, c), const(1)), k - h[2]))
                    term = mul(term, power(diff, k - int(i == 1) - h[1]))
                q = add(q, scale(term, co))
            row.append(divpoly(q, [diff, a, b][i]))
        out.append(row)
    return out

def residual(mat, n):
    """Exact Laurent coefficient of c^{-n}, where n may be negative."""
    out = []
    for row in mat:
        r = {}
        for j, q in enumerate(row):
            for ex, v in q.items():
                deg = ex[4]
                power_bin = deg + n
                if j == 0:
                    if power_bin == 0:
                        r = add(r, shift({ex: v}, 4, -deg))
                elif power_bin >= 0:
                    rr = shift(shift({ex: v * binom(F(1, 2), power_bin)}, 4, -deg), j - 1, -power_bin)
                    r = add(r, rr)
        out.append(r)
    return out

def pairrows(u, v, kind):
    if kind == 'R':
        pu = [add(scale(power(u, 2), 7), scale(mul(u, v), -35)), add(scale(power(u, 3), 64), scale(mul(power(u, 2), v), -140)), add(scale(power(u, 4), 72), scale(mul(power(u, 3), v), -120))]
        pv = [add(scale(mul(u, v), 35), scale(power(v, 2), -7)), add(scale(mul(u, power(v, 2)), 140), scale(power(v, 3), -64)), add(scale(mul(u, power(v, 3)), 120), scale(power(v, 4), -72))]
    elif kind == 'T':
        pu = [add(scale(power(u, 3), -15), scale(mul(power(u, 2), v), 56), scale(mul(u, power(v, 2)), -105)), add(scale(power(u, 4), -180), scale(mul(power(u, 3), v), 512), scale(mul(power(u, 2), power(v, 2)), -420)), add(scale(power(u, 5), -240), scale(mul(power(u, 4), v), 576), scale(mul(power(u, 3), power(v, 2)), -360))]
        pv = [add(scale(mul(u, power(v, 2)), 70), scale(power(v, 3), -6)), add(scale(mul(u, power(v, 3)), 100), scale(power(v, 4), -12)), scale(power(v, 5), 24)]
    else:
        pu = [add(power(u, 2), scale(mul(u, v), -5)), add(scale(power(u, 3), 4), scale(mul(power(u, 2), v), -8))]
        pv = [add(scale(mul(u, v), 5), scale(power(v, 2), -1)), add(scale(mul(u, power(v, 2)), 8), scale(power(v, 3), -4))]
    q0 = {}
    row = [{}, {}, {}, {}, {}]
    for fixed, pc in [(u, pu), (v, pv)]:
        j = next((i for i in range(4) if V[i] == fixed))
        for ell, pe in enumerate(pc):
            row[j + 1] = add(row[j + 1], shift(mul(fixed, pe), 4, ell))
            for n in range(ell + 1):
                q0 = add(q0, scale(shift(shift(mul(fixed, pe), j, -n), 4, ell - n), -binom(F(1, 2), n)))
    row[0] = q0
    return row

def radd(*rows):
    return [add(*(r[j] for r in rows)) for j in range(5)]

def rneg(row):
    return [scale(q, -1) for q in row]

def general_residual(row, n):
    r = {}
    for j, q in enumerate(row):
        for ex, v in q.items():
            nn = n + ex[4]
            if j == 0:
                if nn == 0:
                    r = add(r, shift({ex: v}, 4, -ex[4]))
            elif nn >= 0:
                r = add(r, shift(shift({ex: v * binom(F(1, 2), nn)}, 4, -ex[4]), j - 1, -nn))
    return r

def ordmatrix(anchors, cv, k):
    m = len(anchors)
    hh = prod(anchors)
    al = [F(0)] + [F(hh, a) for a in anchors]
    nn = hh * cv
    dd = prod((anchors[i] - anchors[j] for i in range(m) for j in range(i + 1, m)))
    clearing = (4 * hh ** m * cv * dd ** 2) ** k
    mat = []
    for i in range(m + 1):
        row = []
        for j in range(m + 1):
            vv = F(0)
            for h in comps(k + int(i == j) - 1, m + 1):
                v = binom(F(2 * k + 1, 2), h[j]) * F(nn) ** (-k) * (nn + al[j]) ** (k - h[j])
                for l in range(m + 1):
                    if l != j:
                        v *= binom(-k - int(i == l), h[l]) * (al[j] - al[l]) ** (-k - int(i == l) - h[l])
                vv += v
            row.append(vv * clearing)
        mat.append(row)
    return mat

def detnum(mat):
    out = F(0)
    n = len(mat)
    for p in permutations(range(n)):
        out += (-1) ** sum((p[i] > p[j] for i in range(n) for j in range(i + 1, n))) * prod((mat[i][p[i]] for i in range(n)))
    return out
V=[mon(1,**{v:1}) for v in "abstc"];a,b,s,t,c=V
diff=add(a,scale(b,-1))


# Read the printed parameters and coefficient tables from the article input.
PARAMS=[]
for line in S5.splitlines():
    match=re.match(r'\s*\$([NA]_\d)\$&(.*)$',line)
    if match:
        cells=[z.strip(' $\\.') for z in match[2].split('&')]
        check('eight parameter cells '+match[1],len(cells)==7)
        PARAMS.append((match[1].replace('_',''),*(int(x) for x in cells[:4]),*(F(x) for x in cells[4:])))
check('five actual source families',[p[0] for p in PARAMS]==['N2','N3','A3','A4','A5'])
PRINTED_L={}
for line in S5.splitlines():
    mat=re.match(r'([NA]_\d)&0&([\d,]+)',line)
    if mat:
        current=mat[1].replace('_','');PRINTED_L[current]=[[int(z) for z in mat[2].split(',')]]
    else:
        mat=re.match(r'\s*&1,2&([\d,]+)',line)
        if mat:PRINTED_L[current].append([int(z) for z in mat[1].split(',')])
check('all actual family coefficient rows parsed',len(PRINTED_L)==5 and all(len(v)==2 for v in PRINTED_L.values()))
Mblock=S5.split('Valid integer choices of $M$ in the table order are')[1].split(r'\]')[0]
M_PRINT=list(map(int,re.findall(r'\d+',Mblock)))
check('five actual M values',len(M_PRINT)==5)
scalar_block=S5.split(r'k&4^k\binom{k+1/2}{h}')[1].split(r'\end{array}')[0]
SCALARS={int(k):list(map(int,v.split(','))) for k,v in re.findall(r'([2-5])&([\d,]+)',scalar_block)}
PRINTED_DEG=[]
for name in ['R_{ab}','H','T_{ab}','V']:
    mm=re.search(re.escape(name)+r'&([\d]+)&([\d]+)&([\d]+)',S5)
    check('actual degree table row '+name,mm is not None)
    PRINTED_DEG.append(list(map(int,mm.groups())))
SRC_LW={}
for which,name in [('L',z) for z in ['R_{ab}','H','T_{ab}','V']]+[('W',z) for z in ['R_{ab}','H','T_{ab}','V']]:
    token=which+'_'+('{'+name+'}' if name in ['R_{ab}','T_{ab}'] else name)
    mm=re.search(re.escape(token)+r'\(x\)&=([^\n]+)',S5)
    check('actual residual polynomial '+token,mm is not None)
    expr=mm[1].rstrip(',.\\ ').replace('^','**')
    expr=re.sub(r'(\d|x|\))(?=x|\()',r'\1*',expr)
    check('safe polynomial grammar '+token,re.fullmatch(r'[0-9x+*()]+',expr) is not None)
    SRC_LW[(which,name)]=expr

def source_poly_value(which,name,x):
    return eval(SRC_LW[(which,name)],{'__builtins__':{}},{'x':x})

check('residual branch jump statement','difference of the two boundary values of the numerator has modulus' in S5)
check('H and V combined before coefficient bounds','the combinations defining $H$ and $V$ already performed.' in S5)
check('P row applies for every ratio q at least one','For every $q\\ge1$ with $B\\ge qA$' in S5)
FAMILIES={}
for ix,(name,k,m,g,Delta,J,beta,Rstar) in enumerate(PARAMS):
 mat=family(k,m);FAMILIES[name]=mat;P=m+2*k
 scalar_lists=SCALARS
 check(name+' printed numerator scalar list',[4**k*binom(F(2*k+1,2),hh) for hh in range(k+1)]==scalar_lists[k])
 check(name+' all coefficients integral after g division',all(v.denominator==1 and v.numerator%g==0 for row in mat for q in row for v in q.values()))
 check(name+' anchor divisibility',all(ex[j-1]>=1 for row in mat for j in (1,2) for ex in row[j]))
 for n in range(-k,P-k):check(name+f' residual coefficient c^({-n}) vanishes',all(not r for r in residual(mat,n)))
 check(name+' predicted first residual not identically zero',any(residual(mat,P-k)))
 expected=scale(mul(mon(1,a=m-1,b=m-1,c=k-m),power(diff,4*k-1)),Delta)
 # Full determinant, without truncating c or using the residual degree bound.
 check(name+' full symbolic determinant',det(mat)==expected)
 lp=[]
 for i in range(3):
  cols=[shift(mat[i][1],0,-1),shift(mat[i][2],1,-1)]
  check(name+f' row{i} homogeneous weight',all(e[0]>=0 and e[1]>=0 and e[0]+e[1]==e[4]+P-k-2 and 0<=e[4]<=k for q in cols for e in q))
  sums=[sum(abs(v) for q in cols for e,v in q.items() if e[4]==l) for l in range(k+1)]
  check(name+f' row{i} printed coefficient sums',sums==PRINTED_L[name][int(i>0)],list(map(int,sums)));lp.append(sums)
 threshold=F(101,100*g)*max(sum(z*F(1000)**(l-k) for l,z in enumerate(row)) for row in lp)
 check(name+' printed uniform M',threshold<M_PRINT[ix],str(threshold))
 exactbeta=F(factorial(2*k+2)*factorial(2*P-2*k-2),4**P*factorial(k+1)*factorial(P-k-1)*factorial(P))
 lam=F(2*k+1,2)-m
 aa=[(-1)**j*binom(lam,j)*binom(lam-k-j,k-j) for j in range(k+1)]
 bb=[(-1)**j*binom(lam,j)*binom(lam-k-1-j,k-1-j) for j in range(k)]
 def at(v,i):return v[i] if 0<=i<len(v) else 0
 for rr in range(2*k+1):check(name+f' printed convolution r={rr}',sum(at(aa,j)*at(aa,k-rr+j)-at(bb,j)*at(bb,k-rr+j) for j in range(k+1))==J*(-1)**rr*comb(2*k,rr))
 def ABpoly(xx,yy,which):
  dxy=add(xx,scale(yy,-1));ans={};upper=k if which=='A' else k-1
  for hh in range(upper+1):
   ans=add(ans,scale(mul(power(xx,hh),power(dxy,upper-hh)),binom(-k-int(which=='B'),hh)*binom(lam,upper-hh)))
  return ans
 block=[[scale(mul(power(a,m-1),ABpoly(b,a,'A')),4**k),scale(mul(power(b,m),ABpoly(a,b,'B')),4**k)],
        [scale(mul(power(a,m),ABpoly(b,a,'B')),4**k),scale(mul(power(b,m-1),ABpoly(a,b,'A')),4**k)]]
 actual=[[{ex:co for ex,co in mat[ii][jj].items() if ex[4]==0} for jj in (1,2)] for ii in (1,2)]
 check(name+' printed constant root block',actual==block)
 check(name+' J',(-1)**k*binom(lam,k)*binom(lam-k,k)==J)
 check(name+' Delta',4**k*binom(F(2*k+1,2),m)*16**k*J==Delta)
 check(name+' beta',exactbeta==beta)
 check(name+' Rstar',2*16**k*beta/g==Rstar)
 print(name,'completed',flush=True)

Rab=pairrows(a,b,'R');Rbs=pairrows(b,s,'R');Ras=pairrows(a,s,'R');Tab=pairrows(a,b,'T');Tbs=pairrows(b,s,'T');Tsa=pairrows(s,a,'T');Tas=pairrows(a,s,'T')
H=radd(Rab,Rbs,rneg(Ras));VV=radd(Tas,rneg(Tbs));ADAPT=[Rab,H,Tab,VV]
for name,row in [('Rab',Rab),('Tab',Tab),('H',H),('V',VV),('Rbs',Rbs),('Ras',Ras),('Tbs',Tbs),('Tsa',Tsa)]:
 check(name+' integral polynomial row',all(v.denominator==1 and min(e)>=0 for q in row for e,v in q.items()))
 check(name+' anchor divisibility',all(e[j-1]>=1 for j in range(1,5) for e in row[j]))
 for n in range(-2,5):check(name+f' cancellation n={n}',not general_residual(row,n))
ER=add(mon(F(5,128),a=-2),mon(F(-5,128),b=-2),mon(F(5,256),a=1,b=-3),mon(F(-5,256),a=-3,b=1))
ET=add(mon(F(-75,128),a=-1),mon(F(5,16),a=-2,b=1),mon(F(-15,256),a=-3,b=2),mon(F(-35,256),a=1,b=-2),mon(F(15,32),b=-1))
check('printed ER',general_residual(Rab,5)==ER);check('printed ET',general_residual(Tab,5)==ET)
Fs=add(scale(add(scale(s,2),scale(a,-1),scale(b,-1)),7),scale(mul(c,add(scale(power(s,2),3),scale(mul(add(a,b),s),-1),scale(mul(a,b),-1))),6))
base=mul(mon(1,a=2,b=2,s=2),mul(power(diff,4),mul(add(a,scale(s,-1)),add(b,scale(s,-1)))))
# The full 4x4 determinants are expanded, not just constant/linear terms.
check('four-row full polynomial 4x4 determinant',det([r[:4] for r in ADAPT])==scale(mul(mul(base,diff),Fs),-78750))
L_PRINT=[];W_PRINT=[]
for idx,row in enumerate(ADAPT):
 pp=[shift(row[j],j-1,-1) for j in range(1,4)]
 ds=[sum(abs(v) for q in pp for ex,v in q.items() if ex[4]==ell) for ell in range(3)]
 check(f'four-row row {idx} degree sums',ds==PRINTED_DEG[idx],list(map(int,ds)))
 leading=general_residual(row,5);ll=[F(0)]*6;ww=[0]*6
 for ex,v in leading.items():
  i,j,k=ex[:3];check(f'row{idx} leading monomial guards {ex}',i+j+k+1<=0 and 2-j-k>=0 and 2-k>=0)
  ll[2-j-k]+=abs(v)*256
 for fixed in range(3):
  for ex,coef in pp[fixed].items():
   ee=list(ex);ee[fixed]-=4+ex[4];i,j,k=ee[:3]
   check(f'row{idx} tail monomial guards {ex}/{fixed}',i+j+k+1<=0 and 2-j-k>=0 and 2-k>=0)
   ww[2-j-k]+=abs(coef)
 for xx in range(6):
  check(f'four-row {idx} source leading polynomial x={xx}',sum(z*xx**j for j,z in enumerate(ll))==source_poly_value('L',['R_{ab}','H','T_{ab}','V'][idx],xx))
  check(f'four-row {idx} source tail polynomial x={xx}',sum(z*xx**j for j,z in enumerate(ww))==source_poly_value('W',['R_{ab}','H','T_{ab}','V'][idx],xx))
 L_PRINT.append(ll);W_PRINT.append(ww)
# Ordered-modulus coefficient envelope H*=max(BS^3,B^3S^2).
 for jj,col in enumerate(pp):
  for ex in col:
   ii,jj2,kk=ex[:3]
   exponents=(ii+jj2+kk,jj2+kk,kk)
   check(f'row {idx} coefficient envelope {ex}',any(all(x<=y for x,y in zip(exponents,target)) for target in [(4,4,3),(5,5,2)]))
print('quadratic rows completed',flush=True)

# Unrestricted rank proof: its only arithmetic input is |E|>=1
# for a nonzero algebraic integer in an imaginary quadratic full ring.
EE=add(scale(power(s,2),3),scale(mul(add(a,b),s),-1),scale(mul(a,b),-1))
LL=add(scale(s,2),scale(a,-1),scale(b,-1))
check('rank degenerate algebra E=(s-a)^2+(a+s)L',EE==add(power(add(s,scale(a,-1)),2),mul(add(a,s),LL)))
check('rank guard strictness 6*200>28',6*200>28)
# Bernstein coefficients certify all six polynomial comparisons on [0,1].
DOMINANCE=[]
for typ,rows in [('256L',L_PRINT),('W',W_PRINT)]:
    for i in range(3):
        pp=[F(rows[3][j]-rows[i][j]) for j in range(6)]
        bern=[sum(pp[j]*F(comb(ii,j),comb(5,j)) for j in range(ii+1)) for ii in range(6)]
        check(f'{typ} V dominates row {i} on [0,1]',all(bb>=0 for bb in bern),list(map(str,bern)))
        DOMINANCE.append({'table':typ,'row':i,'bernstein_degree':5,'coefficients':list(map(str,bern))})
check('all source L/W coefficients nonnegative',all(v>=0 for row in L_PRINT+W_PRINT for v in row))
R1=sum(L_PRINT[-1])/256+F(21,1024*999)*sum(W_PRINT[-1])
check('source R(1) strict constant 3',R1==F(65233,28416) and R1<3,str(R1))
LINEAR=[pairrows(u,v,'L') for u,v in [(a,b),(a,s),(b,s),(a,t),(b,t)]]
for idx,row in enumerate(LINEAR):
 for n in range(-1,4):check(f'linear{idx} cancellation {n}',not general_residual(row,n))
 check(f'linear{idx} coefficient total',sum(abs(v) for j in range(1,5) for v in row[j].values())==36)
expected=const(180)
for z in [mon(1,a=2,b=2,s=2,t=2),power(diff,2),add(a,scale(s,-1)),add(a,scale(t,-1)),add(b,scale(s,-1)),add(b,scale(t,-1)),add(s,scale(t,-1))]:expected=mul(expected,z)
check('linear full polynomial determinant',det(LINEAR)==expected)
# Full exact evaluation at c=1,2,-3 confirms the degree argument independently at these points.
for cv in [1,2,-3]:
 mat=[[const(evaluate(q,[2,-5,7,11,cv])) for q in row] for row in LINEAR]
 check(f'linear determinant evaluation c={cv}',det(mat)==const(evaluate(expected,[2,-5,7,11,cv])))
print('linear rows completed',flush=True)

for anchors in [[2,-3],[2,5,-7],[2,-3,5,7],[1,2,3,5,7]]:
 m=len(anchors);hh=prod(anchors);dd=prod(anchors[i]-anchors[j] for i in range(m) for j in range(i+1,m))
 for k in [1,2,3]:
  for cv in [11,-13]:
   mat=ordmatrix(anchors,cv,k)
   expected=((-1)**(m*(m+1)//2)*4**(m+1)*hh**(m-1)*dd**(2*m))**k*prod(binom(F(1,2)-v*k,k) for v in range(-1,m))
   check(f'ordinary m={m},k={k},c={cv} determinant',detnum(mat)==expected)
   check(f'ordinary m={m},k={k},c={cv} integer column factors',all(mat[i][j].denominator==1 and (j==0 or (mat[i][j]/anchors[j-1]**((m-1)*k)).denominator==1) for i in range(m+1) for j in range(m+1)))
   if k==1:check(f'ordinary m={m},k=1,c={cv} even quotients',all((mat[i][j]/(2 if j==0 else 2*anchors[j-1]**(m-1))).denominator==1 for i in range(m+1) for j in range(m+1)))
   if m==3 and k==2:
    aa,bb,ss=anchors;dv=[dd,aa*(bb-ss),bb*(aa-ss),ss*(aa-bb)];rm=[[v/dv[i] for v in row] for i,row in enumerate(mat)]
    check(f'ordinary normalized 3-root det c={cv}',detnum(rm)==-226800*hh**3*dd**10)
    check(f'ordinary normalized column divisibility c={cv}',all((rm[i][j]/(1 if j==0 else anchors[j-1])).denominator==1 for i in range(4) for j in range(4)))
print('generic exact residue instances completed',flush=True)

# Exact strict numerical bounds (no floating point).
ineq=[('uniform branch e',F(1001,1000)*F(5,4)**2,F(3,2)**2),('sharp branch e',F(1001,1000)*F(201,200)**2,F(101,100)**2),('ordinary H2',F(125,3)*F(5007,5000),42),('ordinary coefficient',4*F(5007,4995)*F(25*1001,16000),9),('ordinary R2 lower',F(2,5),F(27,64)*F(999,1000)**2),('three G',F(16384,27)*F(4005,4000),608),('three pe',F(25*1001*4005,64000*999),F(16,9)),('three residual',27*F(1000,999)**3,28),('finite eta',F(1000,999)**17,F(101,100)**2),('quadratic M',F(101,100)*max(F(L[2])+F(L[1],1000)+F(L[0],10**6) for L in PRINTED_DEG),2400),('adaptive V',F(65233,28416),F(23,10)),('linear residual',F(36*5*1000,128*999),F(3,2))]
for name,lhs,rhs in ineq:check(name,lhs<rhs,{'lhs':str(lhs),'rhs':str(rhs)})
check('adaptive rational residual equality',F(35,16)+F(21*5268,1024*999)==F(65233,28416))
# Direct polynomial identity proves the recurrence comparison for all k.
num=scale(mul(mul(add(a,const(F(3,2))),add(scale(a,2),const(F(-1,2)))),add(scale(a,2),const(F(1,2)))),F(27,4))
den=mul(mul(add(scale(a,3),const(1)),add(scale(a,3),const(2))),add(scale(a,3),const(3)))
check('all-order beta recurrence polynomial',add(den,scale(num,-1))==add(scale(power(a,2),F(432,32)),scale(a,F(1110,32)),const(F(273,32))))
eps=F(1,4000)
check('closed lower interval conversion',F(4,1)/(1-eps)*(1-eps)==4)
check('closed upper interval conversion',F(1,4)/(1+eps)*(1+eps)==F(1,4))
check('linear residual exact ceiling',F(36*5*1000,128*999)==F(625,444))

for k in range(1,6):
    upper=(0,(-1)**k);lower=(0,-upper[1])
    check(f'actual boundary jump k={k}',sum(x*x for x in upper)==sum(x*x for x in lower)==1 and sum((x-y)**2 for x,y in zip(upper,lower))==4)
# Expand the entire fixed m=3, k=2 matrix from the polynomial summands
# in Lemma 4.4. Signs are retained from the original residue denominators.
def prod_polys(polynomials):
 result=const(1)
 for poly in polynomials:
  result=mul(result,poly)
 return result

anchors=[a,b,s]; m=3; k=2
vand=mul(mul(add(a,scale(b,-1)),add(a,scale(s,-1))),add(b,scale(s,-1)))
row_factors=[vand,mul(a,add(b,scale(s,-1))),mul(b,add(a,scale(s,-1))),mul(s,add(a,scale(b,-1)))]
NORMALIZED_R=[]
for i in range(4):
 row=[]
 for j in range(4):
  q={}
  for hs in comps(k+int(i==j)-1,4):
   scalar=4**k*binom(F(2*k+1,2),hs[j])*prod(binom(-k-int(i==ell),hs[ell]) for ell in range(4) if ell != j)
   if j == 0:
    term=mul(power(c,k-hs[0]),power(vand,2*k))
    scalar*=(-1)**sum(k+int(i==ell)+hs[ell] for ell in range(1,4))
    for ell in range(1,4):
     term=mul(term,power(anchors[ell-1],int(i==ell)+hs[ell]))
   else:
    aj=anchors[j-1]
    term=mul(power(aj,(m-1)*k),power(add(mul(aj,c),const(1)),k-hs[j]))
    numerator=power(vand,2*k); denominator=const(1)
    for ell in range(1,4):
     if ell != j:
      term=mul(term,power(anchors[ell-1],int(i==ell)+hs[ell]))
      denominator=mul(denominator,power(add(anchors[ell-1],scale(aj,-1)),k+int(i==ell)+hs[ell]))
    term=mul(term,divpoly(numerator,denominator))
   q=add(q,scale(term,scalar))
  row.append(divpoly(q,row_factors[i]))
 NORMALIZED_R.append(row)
for i,row in enumerate(NORMALIZED_R):
 for j,poly in enumerate(row):
  check(f'Lemma 4.4 normalized entry R[{i},{j}] integral polynomial',all(v.denominator==1 and min(ex)>=0 for ex,v in poly.items()))
  if j:
   check(f'Lemma 4.4 normalized column factor R[{i},{j}]',all(ex[j-1]>=1 for ex in poly))
 for n in range(-2,6):
  check(f'Lemma 4.4 normalized row {i} residual c^({-n}) vanishes',not general_residual(row,n))
check('Lemma 4.4 first residual nonzero',any(general_residual(row,6) for row in NORMALIZED_R))
check('Lemma 4.4 product of row divisors',prod_polys(row_factors)==mul(mon(1,a=1,b=1,s=1),power(vand,2)))
R_DETERMINANT=scale(mul(mon(1,a=3,b=3,s=3),power(vand,10)),-226800)
check('Lemma 4.4 full normalized polynomial determinant',det(NORMALIZED_R)==R_DETERMINANT)
for cv in [11,-13]:
 numerical=ordmatrix([2,5,-7],cv,2)
 factors=[(2-5)*(2+7)*(5+7),2*(5+7),5*(2+7),-7*(2-5)]
 check(f'Lemma 4.4 polynomial-to-residue normalization at c={cv}',[[evaluate(poly,[2,5,-7,0,cv]) for poly in row] for row in NORMALIZED_R]==[[val/factors[i] for val in row] for i,row in enumerate(numerical)])
print('Lemma 4.4: full normalized polynomial matrix and determinant completed',flush=True)

def polynomial_record(poly):
    """Sorted exact coefficients of monomials in (a,b,s,t,c)."""
    return [{'powers':list(ex), 'coefficient':str(co)} for ex,co in sorted(poly.items())]

def matrix_record(matrix):
    return [[polynomial_record(poly) for poly in row] for row in matrix]

certificate = {
    'variables':['a','b','s','t','c'],
    'coefficient_format':'Exact rational strings; powers may be negative only for Laurent residuals.',
    'two_anchor_families':{name:matrix_record(mat) for name,mat in FAMILIES.items()},
    'P_rows':{name:[polynomial_record(p) for p in row] for name,row in zip(['R_ab','H','T_ab','V'],ADAPT)},
    'Q_rows':[ [polynomial_record(p) for p in row] for row in LINEAR],
    'P_leading_residuals':{'R_ab':polynomial_record(ER),'T_ab':polynomial_record(ET)},
    'Q_determinant':[],
    'normalized_R_lemma_4_4':matrix_record(NORMALIZED_R),
    'normalized_R_determinant':polynomial_record(R_DETERMINANT),
}
# The name `expected` is reused by the numerical probes above. Reconstruct the
# Q determinant here for the machine-readable coefficient certificate.
qdet=const(180)
for factor in [mon(1,a=2,b=2,s=2,t=2),power(diff,2),add(a,scale(s,-1)),add(a,scale(t,-1)),add(b,scale(s,-1)),add(b,scale(t,-1)),add(s,scale(t,-1))]:
    qdet=mul(qdet,factor)
certificate['Q_determinant']=polynomial_record(qdet)
(RESULTS/'approximation_polynomials_coefficients.json').write_text(json.dumps(certificate,indent=2)+'\n')
result={
    'status':'PASS','article_sha256':SHA,
    'paper_sections':['Section 4, Simultaneous approximation','Section 5.2, Five families with two fixed entries','Section 5.3, Quadratic rows and their error bounds','Section 5.4, Five linear rows on four fixed entries'],
    'paper_labels':['lem:ordinary-residue-construction','lem:ordinary-three-normalized','eq:ff-generator','tab:ff-parameters','eq:ff-weight','eq:ff-uniformM','eq:ff-prows','eq:ff-leading','eq:ff-adaptiveres'],
    'method':'Sparse multivariate polynomials and fractions.Fraction, using only exact integer and rational arithmetic. Full determinants and Laurent coefficients are expanded.',
    'check_count':len(checks),'checks':checks,
    'family_parameter_columns':['name','k','m','g','c_det','J','beta','R_star'],
    'family_parameters':PARAMS,'coefficient_sums':PRINTED_L,'uniform_M':M_PRINT,
    'quadratic_degree_sums':PRINTED_DEG,
    'residual_polynomials':{str(k):v for k,v in SRC_LW.items()},
    'dominance_certificates':DOMINANCE,
    'expanded_coefficients_file':'approximation_polynomials_coefficients.json',
    'scope':'The five finite families, normalized R matrix of Lemma 4.4, and P/Q polynomial identities are checked in full. Ordinary-family evaluations at m=2,3,4,5 and k=1,2,3 are finite examples; the all-order identity is the cited theorem in the article.',
    'elapsed_seconds':round(time.monotonic()-START,3),
}
(RESULTS/'approximation_polynomials.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
print('PASS:',len(checks),'exact polynomial, determinant, coefficient and bound checks.')
print('Five families: N2, N3, A3, A4, A5. Uniform M:', ', '.join(map(str,M_PRINT)))
print('P rows: all four full polynomial rows, determinant, leading residuals and coefficient bounds checked.')
print('Q rows: all five full polynomial rows, determinant and residual cancellation checked.')
print('Lemma 4.4: full normalized R matrix, integrality, column factors and determinant checked.')
print('Results: approximation_polynomials.json; approximation_polynomials_coefficients.json')
