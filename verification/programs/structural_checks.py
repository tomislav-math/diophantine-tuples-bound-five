#!/usr/bin/env python3
"""Exact structural identities, the regular-tail enumeration, and gap constants.
All arithmetic is integer, rational, or sparse Laurent-polynomial arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import isqrt
import hashlib, json, re, os, tempfile

# Sparse exact Laurent-polynomial arithmetic. Every division below is by a monomial.
class LP:
    def __init__(self,v=0):
        if isinstance(v,LP):self.d=v.d.copy()
        elif isinstance(v,dict):self.d={tuple(k):F(c) for k,c in v.items() if c}
        else:self.d={():F(v)} if v else {}
    @staticmethod
    def var(name):return LP({((name,1),):1})
    def __hash__(self):return hash(tuple(sorted(self.d.items())))
    def __add__(self,o):
        o=LP(o);d=self.d.copy()
        for m,c in o.d.items():d[m]=d.get(m,0)+c
        return LP(d)
    __radd__=__add__
    def __neg__(self):return LP({m:-c for m,c in self.d.items()})
    def __sub__(self,o):return self+-LP(o)
    def __rsub__(self,o):return LP(o)+-self
    def __mul__(self,o):
        o=LP(o);out={}
        for m,c in self.d.items():
            for n,e in o.d.items():
                z=dict(m)
                for v,k in n:z[v]=z.get(v,0)+k
                key=tuple(sorted((v,k) for v,k in z.items() if k));out[key]=out.get(key,0)+c*e
        return LP(out)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:
            assert len(self.d)==1
            (m,c),=self.d.items();return LP({tuple((v,k*n) for v,k in m):c**n})
        out=LP(1)
        for _ in range(n):out=out*self
        return out
    def __truediv__(self,o):return self*LP(o)**-1
    def __rtruediv__(self,o):return LP(o)*self**-1
    def __eq__(self,o):return self.d==LP(o).d
    def __repr__(self):return str(self.d)
    def coeff(self,v,n):
        name=next(iter(v.d))[0][0];d={}
        for m,c in self.d.items():
            if dict(m).get(name,0)==n:
                k=tuple((w,e) for w,e in m if w!=name);d[k]=d.get(k,0)+c
        return LP(d)
    def subs(self,mapping):
        byname={next(iter(k.d))[0][0]:v for k,v in mapping.items()};out=LP()
        for m,c in self.d.items():
            val=LP(c)
            for v,e in m:val*=byname.get(v,LP.var(v))**e
            out+=val
        return out
class Sym:
    @staticmethod
    def symbols(names,**kwargs):
        z=[LP.var(v) for v in names.split()];return z[0] if len(z)==1 else tuple(z)
    @staticmethod
    def expand(v):return LP(v)
    factor=expand
    Rational=F
    @staticmethod
    def prod(values):
        out=1
        for v in values:out*=v
        return out
    @staticmethod
    def binomial(v,n):
        out=F(1)
        for j in range(n):out*=F(v-j,j+1)
        return out
    class Poly:
        def __init__(self,p,*vs):self.p=p;self.names=[next(iter(v.d))[0][0] for v in vs]
        def terms(self):return [(tuple(dict(m).get(v,0) for v in self.names),c) for m,c in self.p.d.items()]
s=Sym()


ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data/article.txt'
EXPECTED_SOURCE_SHA='9cb6ed2154056e158d2c420012e94ef7c397c69d66c0fd9301a1bb52197dc432'
assert hashlib.sha256(SRC.read_bytes()).hexdigest()==EXPECTED_SOURCE_SHA
if os.environ.get('D1_RESULTS_DIR'):
 OUT=Path(os.environ['D1_RESULTS_DIR']);OUT.mkdir(parents=True,exist_ok=True)
else:
 (ROOT/'reruns').mkdir(exist_ok=True)
 OUT=Path(tempfile.mkdtemp(prefix='enumeration-',dir=ROOT/'reruns'))
source=SRC.read_text();digest=EXPECTED_SOURCE_SHA
checks=[]; results={}
def ck(name,condition,evidence=None):
    assert condition,name
    checks.append({'name':name,'passed':True,'evidence':str(evidence) if evidence is not None else None})

def nz(f):return s.expand(f)==0
# Bind the calculations to the corresponding statements in the article.
for fragment in [r'\label{eq:irregular-quadratic}',r'|1+2ab|-|a+b|',r'If \(A\ge3\), then \(T>3ABS\)',r'F(2)<7/10',r'F(3)<19/20',r'K=400',r'16405041495117',r'\mu=\frac{51051}{400000}']:
    ck('source contains '+fragment,fragment in source)
# Completion identity and the exact quadratic before majorizing coefficients.
a,b,c=s.symbols('a b c'); h=1+2*a*b; q=a+b; delta=(a-b)**2-4
u0=a+b+c+2*a*b*c
uv=s.expand(u0**2-4*(a*b+1)*(a*c+1)*(b*c+1))
ck('regular-completion product',nz(uv-(c*c-2*q*c+delta)),uv)
A,B,X=s.symbols('A B X',nonnegative=True)
fmajor=X*X-2*(2*A*B-1-(A+B))*X+(A+B)**2+4+2*(A+B)+1
coarse=(A+B+X+1)**2+4-4*A*B*X
ck('quadratic majorant identity',nz(fmajor-coarse))
ck('irregular lower endpoint',F(10,3)**2+F(4,9)-12==-F(4,9))
ck('irregular upper endpoint',-1+F(7,9)+F(7,9)**2/4+F(1,81)==-F(19,324))
# Tail analytic constants: branch rationalization inequalities and convex endpoints.
for rho,k in [(F(1,12),F(20,39)),(F(1,3),F(5,9)),(F(1,4),F(7,13)),(F(9,25),F(5,9))]:
    ck('branch rationalization '+str(rho),(1/k-1)**2<=1-rho)
q0=F(7,13)+F(125,216); r0=F(625,1053)
ck('tail endpoint at A=1',F(35,108)+q0/2+F(7,48)*r0==F(12247,12636)<1)
ck('tail endpoint at A=B',F(5,27)+F(7,24)*q0+r0/4==F(14813,22464)<1)
ck('tail relative loss',F(35,468)<F(1,13))
ck('tail p lower bound',F(12,13)>F(24,25)**2)
ck('tail perturbation radius',(1+F(35,117))/F(48,13)==F(19,54)<F(9,25))
# Fresh finite enumeration, exactly the domain described in the final tail lemma.
def squarefree(n):return all(n%(p*p) for p in range(2,isqrt(n)+1))
fields=[d for d in range(1,65) if squarefree(d) and (d if d%4==3 else 4*d)<=64]
ck('tail field list',fields==[1,2,3,5,6,7,10,11,13,14,15,19,23,31,35,39,43,47,51,55,59])
def params(d):return (1,(d+1)//4) if d%4==3 else (0,d)
def ring(d):
    tau,nu=params(d)
    def norm(z):x,y=z;return x*x+tau*x*y+nu*y*y
    def add(z,w):return(z[0]+w[0],z[1]+w[1])
    def mul(z,w):x,y=z;u,v=w;return(x*u-nu*y*v,x*v+y*u+tau*y*v)
    def neg(z):return(-z[0],-z[1])
    def conj(z):return(z[0]+tau*z[1],-z[1])
    def disk(n):
        out=[];disc=4*nu-tau*tau
        for y in range(-isqrt(4*n//disc),isqrt(4*n//disc)+1):
            bound=isqrt(4*n-disc*y*y)
            for z in range(-bound,bound+1):
                if (z-tau*y)%2==0:out.append(((z-tau*y)//2,y))
        return sorted(out,key=lambda z:(norm(z),z))
    return norm,add,mul,neg,conj,disk
counts=[0,0,0,0];candidate_records=[];orbits=set();perfield=[]
for d in fields:
    norm,add,mul,neg,conj,disk=ring(d);small=disk(15);roots=disk(70);squares={}
    for r in roots:squares.setdefault(mul(r,r),[]).append(r)
    bc=[z for z in small if 9<=norm(z)<=15];local=[0,0,0,0]
    for b,c in combinations(bc,2):
        prod=add(mul(b,c),(1,0))
        if prod not in squares:continue
        local[0]+=1
        for r in squares[prod]:
            t=add(add(b,c),(2*r[0],2*r[1]))
            if t in [(0,0),b,c] or norm(t)<norm(c):continue
            local[1]+=1
            for a in small:
                if a in [(0,0),b,c,t] or norm(a)>norm(b):continue
                if add(mul(a,b),(1,0)) not in squares or add(mul(a,c),(1,0)) not in squares:continue
                local[2]+=1;val=add(mul(a,t),(1,0));passed=val in squares
                local[3]+=int(passed)
                candidate_records.append({'d':d,'a':a,'b':b,'c':c,'t':t,'at_plus_1':val,'final_pass':passed})
                orbit=[]
                for sign,conjugate,swap in product([False,True],repeat=3):
                    entries=[a,b,c,t]
                    if sign:entries=[neg(z) for z in entries]
                    if conjugate:entries=[conj(z) for z in entries]
                    if swap:entries[1],entries[2]=entries[2],entries[1]
                    orbit.append(tuple(entries))
                orbits.add((d,min(orbit)))
    counts=[x+y for x,y in zip(counts,local)];perfield.append({'d':d,'counts':local})
ck('tail enumeration counts',counts==[90,118,16,0],counts)
ck('tail candidate orbit count',len(orbits)==6,len(orbits))
results['tail_enumeration']={'counts':counts,'per_field':perfield,'candidates':candidate_records,'canonical_orbits':sorted(orbits)}

printed=[]
for line in source.splitlines():
    if re.fullmatch(r'(2|3|6|11|35)&\([^\n]+',line):
        cells=line.rstrip('\\').split('&')
        if len(cells)==6 and all(re.fullmatch(r'\(-?\d+,-?\d+\)',x) for x in cells[1:]):
            printed.append((int(cells[0]),*[tuple(map(int,x[1:-1].split(','))) for x in cells[1:]]))
printed_orbits=set()
for d,a,b,c,t,val in printed:
    norm,add,mul,neg,conj,disk=ring(d)
    ck('printed tail final value '+str((d,a,b,c,t)),add(mul(a,t),(1,0))==val)
    variants=[]
    for sign,conjugate,swap in product([False,True],repeat=3):
        entries=[a,b,c,t]
        if sign:entries=[neg(z) for z in entries]
        if conjugate:entries=[conj(z) for z in entries]
        if swap:entries[1],entries[2]=entries[2],entries[1]
        variants.append(tuple(entries))
    printed_orbits.add((d,min(variants)))
ck('printed six tail rows exactly cover recomputed orbits',len(printed)==6 and printed_orbits==orbits)

# Small determinant estimates and rank-two algebra.
e=F(128,8649);sm=F(1,2)+3*e;hm=F(3,2)+3*e
ck('product-four first constant',hm*sm**2*(1+F(3,36)+sm/324)<F(1,2))
ck('product-four weighted maximum',max(F(1,9),F(11,96))<F(1,8))
s2=F(59,96);h2=F(155,96)
ck('product-four second constant',h2*s2**2*(F(13,12)+s2/324)==F(18212679025,27518828544)<F(2,3))
ck('product-five inverse residual',F(2,3)+F(8,3)*(F(1,200)+F(1,12))==F(203,225)<F(11,12))
ck('product-five determinant prebound',F(3,16)+F(3,2)*F(7,192)+9*F(7,192)**2<F(33,128))
ck('product-five final bound',F(33,128)*(1+F(11,54)+F(505,11664)+F(1,2916))<F(1,3))
ck('product-five nonsingularity bound',F(1,8)+F(67,300)+F(201,200)**2/2916<F(1,2))
h,p0,q0,d=s.symbols('h p0 q0 d')
rank2=(1-h/d)*(1-h/d-p0/d**2)-(-p0/d)*(-q0/d-h/d**2)
ck('rank-two determinant identity',s.factor(rank2-(1-2*h/d+(h*h-p0*q0-p0)/d**2))==0)
# Independent 24-permutation expansion, not a stored polynomial or old check.
edges=list(combinations(range(4),2));zz=s.symbols('e12 e13 e14 e23 e24 e34');edge=dict(zip(edges,zz))
def ent(i,j):return 1 if i==j else 1+edge[tuple(sorted((i,j)))]
def parity(perm):return (-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
det=s.expand(sum(parity(pp)*s.prod(ent(i,pp[i]) for i in range(4)) for pp in permutations(range(4))))
poly=s.Poly(det,*zz);cubic=sum(co*s.prod(z**e for z,e in zip(zz,mon)) for mon,co in poly.terms() if sum(mon)==3)
quartic=sum(co*s.prod(z**e for z,e in zip(zz,mon)) for mon,co in poly.terms() if sum(mon)==4)
triangles=sum(s.prod(edge[tuple(sorted(x))] for x in combinations(t,2)) for t in combinations(range(4),3))
paths=set(tuple(sorted(tuple(sorted((pp[j],pp[j+1]))) for j in range(3))) for pp in permutations(range(4)))
pathsum=sum(s.prod(edge[x] for x in path) for path in paths)
squareedge=0
for i,j in edges:
    k,l=[v for v in range(4) if v not in (i,j)];squareedge+=edge[(i,j)]**2*edge[(k,l)]
ck('cubic determinant graph formula',nz(cubic-(2*triangles-2*pathsum+2*squareedge)))
matching=[edge[(0,1)]*edge[(2,3)],edge[(0,2)]*edge[(1,3)],edge[(0,3)]*edge[(1,2)]]
ck('quartic determinant formula',nz(quartic-(sum(v*v for v in matching)-2*sum(v*w for v,w in combinations(matching,2)))))
ck('no other determinant degrees',all(sum(mon) in [3,4] for mon,co in poly.terms()))
anchors=s.symbols('a b s t',nonzero=True);sub={edge[(i,j)]:1/(2*anchors[i]*anchors[j]) for i,j in edges};prodanchors=s.prod(anchors)
lead=s.factor(prodanchors*cubic.subs(sub));target=sum(anchors[i]**2/prodanchors for i in range(4))/4-sum(1/(anchors[i]*anchors[j]) for i,j in edges)/4
ck('substituted leading Laurent identity',s.factor(lead-target)==0)

grouped={}
for mon,co in s.Poly(cubic,*zz).terms():
    powers=tuple(1-sum(mon[j] for j,e in enumerate(edges) if i in e) for i in range(4))
    grouped[powers]=grouped.get(powers,F(0))+abs(co)
ck('cubic absolute coefficient groups',len(grouped)==10 and sorted(grouped.values())==[2]*4+[6]*6)
ck('quartic absolute coefficient sum',sum(abs(co) for mon,co in s.Poly(quartic,*zz).terms())==9)
results['determinant']={'cubic_monomials':len(s.Poly(cubic,*zz).terms()),'quartic_monomials':len(s.Poly(quartic,*zz).terms()),'grouped_absolute_coefficients':{str(k):str(v) for k,v in grouped.items()}}

for mon,co in s.Poly(cubic,*zz).terms():ck('cubic incident count '+str(mon),sum(mon[j] for j,e in enumerate(edges) if 0 in e)<=2)
for mon,co in s.Poly(quartic,*zz).terms():ck('quartic incident count '+str(mon),sum(mon[j] for j,e in enumerate(edges) if 0 in e)==2)
rho=s.symbols('rho');ck('small defect polynomial',nz((1-rho)*(7+2*rho)**2-(7-2*rho)**2-rho*(7-28*rho-4*rho*rho)))
ck('small defect endpoints',7-28*F(1,5)-4*F(1,5)**2>0 and (F(1,2)+F(1,35))**2<F(2,7))
def majorant(alpha,beta,tau):
    k1=F(1,2)+1/(7*alpha*beta);k2=F(1,2)+1/(7*beta*beta);K=k1*k1*k2;L=6*K-F(1,2);Q=9*k1*k1*k2*k2
    P=2/(alpha*beta)+1/(beta*beta)
    J=2*K*(2/(alpha*alpha*beta*beta)+1/beta**4)+L*(1/(alpha*alpha*beta*beta)+2/(alpha*beta**3))+Q/(alpha*alpha*beta**4)
    ck('majorant convexity '+str((alpha,beta,tau)),J>0)
    return 2*K*tau+L*P+J/tau
rows=[(F(1),F(8),F(27,8)),(F(140,99),F(4),F(3)),(F(12,7),F(3),F(11,4)),(F(2),F(3),F(29,10)),(F(20,9),F(3),F(3))]
front=[]
for alpha,beta,t in rows:
    lo=majorant(alpha,beta,F(2,3));hi=majorant(alpha,beta,t)
    ck('frontier AB guard '+str(alpha),alpha*beta>5)
    ck('frontier endpoint '+str(alpha),lo<F(1,2) and hi<F(499,500))
    front.append({'alpha':str(alpha),'B0':str(beta),'t':str(t),'F2_3':str(lo),'Ft':str(hi)})
ck('frontier maximum',max(F(row['Ft']) for row in front)==F(964921210337222547,967173115740160000))
for val in [2,3]:
    out=majorant(F(3),F(3),F(val));ck('all-large F('+str(val)+')',out<({2:F(7,10),3:F(19,20)}[val]),out)
results['frontier']=front

# A cleared rational-function identity verifies the barycentric formula for all n,
# rather than checking a finite list of Taylor coefficients.
W=s.symbols('W')
for m in (2,3,4):
    aa=s.symbols(' '.join('a'+str(j) for j in range(m)))
    vand=s.prod(aa[i]-aa[j] for i,j in combinations(range(m),2))
    den=s.prod(1-W/v for v in aa)
    lhs=0
    for i in range(m):
        quotient=(-1)**i*s.prod(aa[j]-aa[k] for j,k in combinations(range(m),2) if i not in (j,k))
        lhs+=quotient*aa[i]**(m-1)*s.prod(1-W/aa[j] for j in range(m) if j!=i)
    rhs=vand*den+(-1)**(m-1)*vand*W**m/s.prod(aa)
    ck('all-order barycentric generating identity m'+str(m),nz(lhs-rhs))
a,b,p,q=s.symbols('a b p q')
relation=(a+b-p-q)**2-4*(a*b+1)*(p*q+1)
eliminated=(2*b*p*q+p+q-a+b)**2-4*(p*q+1)*(b*p+1)*(b*q+1)
ck('strong-two squared elimination identity',nz(relation-eliminated))
mid=a+b+p+2*a*b*p
ck('strong-two quadratic discriminant identity',nz(mid**2-((a+b-p)**2-4*(a*b+1))-4*(a*b+1)*(a*p+1)*(b*p+1)))

# Interpolation low-degree coefficients are multiplied afresh.
x,y,W=s.symbols('x y W');series=sum(s.binomial(s.Rational(1,2),j)*(x*W)**j for j in range(5))*sum(s.binomial(s.Rational(1,2),j)*(y*W)**j for j in range(5))
forms={2:-(x-y)**2/8,3:(x*x-y*y)*(x-y)/16,4:-(x-y)**2*(5*x*x+6*x*y+5*y*y)/128}
for n,target in forms.items():ck('interpolation coefficient '+str(n),nz(s.expand(series).coeff(W,n)-target))
parameters=[]
for H,lam,U in [(32,F(1,10),F(179,10)),(48,F(19,200),F(15))]:
    r=F(5,32)*((F(H,H-1))**3-1);lo=F(45,1024)-r;hi=F(85,1024)+r
    ck('two-anchor remainder H'+str(H),F(200,99)*(F(25,128)+F(1,4)*(F(H,H-1)**2-1))<F(1,2))
    ck('three-anchor residual H'+str(H),lo>0 and F(100,99)*hi<lam)
    ck('three-anchor vanishing window H'+str(H),F(100,201)/hi>5 and F(100,199)/lo<U)
    parameters.append({'H':H,'lower_Q':str(F(100,201)/hi),'upper_Q':str(F(100,199)/lo)})
ck('two-anchor main coefficient',F(202,199)<F(21,20))
ck('three-anchor main coefficient',F(101,199)<F(51,100))
ck('strong-two zero exclusion',11776<13456)
ck('strong-two first term',F(21,20)**2/F(159)<F(1,4))
ck('general-three coefficients',F(296,71)<5 and F(6561,2240)<3)
r=F(5,32)*(F(1000,999)**3-1);lo=F(45,1024)-r;hi=F(85,1024)+r
ck('separated-three estimates',F(808,199)<5 and F(800,99)*hi<F(3,4) and F(100,201)/hi>5 and F(100,199)/lo<15)
ck('four-anchor nonvanishing gap',32-1>F(121,4))
k=F(1,1000);r=F(7,64)*(F(1000,999)**4-1);lo=(1-k)**2*(5-6*k-5*k*k)/128-r;hi=(1+k)**2*(5+6*k+5*k*k)/128+r
ck('direct-four estimates',lo>0 and F(2424,199)<13 and F(2400,99)*hi<F(49,50) and F(100,201)/hi>10 and F(100,199)/lo<17)
results['interpolation_parameters']=parameters
# Retained gap inequalities, retaining strictness at all minimal hypotheses.
h0=F(10001,10000);k=F(3999,4000);mu=F(51051,400000)
ck('retained nonvanishing window',F(179,10)*h0<18)
ck('retained W first contribution',(F(5,4)**2*h0/k)<F(5005,4000)**2)
ck('retained mu',F(51,500)*F(5005,4000)==mu)
alpha_min=10*(F(1,8)-mu/9)
ck('retained alpha and loss',alpha_min>1 and F(9999,10000)*alpha_min**2>1)
ck('retained direct-four first contribution',F(169)/(16**6*k*3**4)<F(1,10000))
large=3**5*F(9999,10000)*F(200,19)**2*(F(1,8)-mu/81)**2
ck('K400 exact value',large==F(16405041495117,40000000000)>400,large)
results.update({'source':'data/article.txt','source_sha256':digest,'assertions':len(checks),'checks':checks,'outcome':'PASS'})
(OUT/'structural_checks.json').write_text(json.dumps(results,indent=2)+'\n')
print('PASS',len(checks),'exact assertions')
print('Regular-tail enumeration: 90 pairs,118 triples,16 candidates,6 orbits,0 extensions')
for row in front:print('Phase-bound row',json.dumps(row,sort_keys=True))
print('Article SHA256',digest)
print('Detailed results: structural_checks.json')
