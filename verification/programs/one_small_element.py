#!/usr/bin/env python3
"""Exact checks for the one-small-entry case in Sections 8 and 9.

Regenerates the finite-family coefficients, enumerates 40 pairs and 57 triples,
checks all five regions and all 41 fixed-pair certificates, and reconstructs
all 56 triple cutoffs plus the separate (1,3,8) cutoff. No result cache is read.
Python standard library only. Floating-point diagnostics are display only.
"""
import sys
if sys.flags.optimize:
    raise SystemExit("Run without -O or -OO: the verification requires assertions.")
from fractions import Fraction as F
from math import isqrt,comb
from pathlib import Path
from itertools import permutations
import json,re,hashlib,os

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/article.tex'
TEX=SOURCE.read_text(encoding='utf-8')
SOURCE_HASH=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
EXPECTED='3653744c838aa8daae33a6e0ce52bf64034f14142af27a2343af7fa6610fe95f'
assert SOURCE_HASH==EXPECTED
import tempfile
if os.environ.get('D1_RESULTS_DIR'):
    OUT = Path(os.environ['D1_RESULTS_DIR'])
    OUT.mkdir(parents=True, exist_ok=True)
else:
    (ROOT/'reruns').mkdir(exist_ok=True)
    OUT = Path(tempfile.mkdtemp(prefix=Path(__file__).stem+'-', dir=ROOT/'reruns'))
    print('Results directory:', OUT)
COUNTER=0
def require(x,msg):
    global COUNTER
    COUNTER+=1
    if not x: raise AssertionError(msg)
def floorroot(v,j=0):
    v=F(v); L=10**j
    return F(isqrt((v.numerator*L*L)//v.denominator),L)
def ceilroot(v,j=0):
    a=floorroot(v,j)
    return a if a*a==v else a+F(1,10**j)

class Ring:
    def __init__(self,d):
        self.d=d; self.t=int(d%4==3); self.u=(d+1)//4 if self.t else d
        self.D=4*self.u-self.t*self.t
    def add(self,a,b): return a[0]+b[0],a[1]+b[1]
    def sub(self,a,b): return a[0]-b[0],a[1]-b[1]
    def neg(self,a): return -a[0],-a[1]
    def mul(self,a,b):
        x,y=a;v,w=b
        return x*v-self.u*y*w,x*w+y*v+self.t*y*w
    def power(self,a,n):
        o=(1,0)
        for _ in range(n): o=self.mul(o,a)
        return o
    def norm(self,a):
        x,y=a
        return x*x+self.t*x*y+self.u*y*y
    def conj(self,a): return a[0]+self.t*a[1],-a[1]
    def ball(self,L):
        yl=isqrt(4*L//self.D)
        for y in range(-yl,yl+1):
            el=isqrt(4*L-self.D*y*y)
            for e in range(-el,el+1):
                if (e-self.t*y)%2==0:
                    yield ((e-self.t*y)//2,y)
    def square(self,z):
        nz=self.norm(z); n=isqrt(nz)
        if n*n!=nz:return False
        p=2*z[0]+self.t*z[1];q=z[1]
        e2=2*n+p; f2num=2*n-p
        if e2<0 or f2num<0 or f2num%self.D:return False
        e=isqrt(e2);f=isqrt(f2num//self.D)
        if e*e!=e2 or self.D*f*f!=f2num:return False
        if e and e*f!=q: f=-f
        if e*f!=q or (e-self.t*f)%2:return False
        rt=((e-self.t*f)//2,f)
        require(self.mul(rt,rt)==z,'square witness failed')
        return True
    def orbitkey(self,a,b):
        aa,ab=a;ba,bb=b
        return (2*aa+self.t*ab<0,ab<0,2*ba+self.t*bb<0,bb<0,aa,ab,ba,bb)
    def canonical(self,a,b):
        p=[(a,b),(self.neg(a),self.neg(b)),(self.conj(a),self.conj(b)),(self.neg(self.conj(a)),self.neg(self.conj(b)))]
        return min(p,key=lambda z:self.orbitkey(*z))

pairsec=TEX.split('The forty first-pair representatives')[1].split('\\end{longtable}')[0]
pairs={}
for m in re.finditer(r'^(\d+)&(\d+)&\$\((-?\d+),(-?\d+)\)\$&\$\((-?\d+),(-?\d+)\)\$&(\d+)&(\d+)',pairsec,re.M):
    i,d,ax,ay,bx,by,na,nb=map(int,m.groups());pairs[i]=(d,(ax,ay),(bx,by))
    require(Ring(d).norm((ax,ay))==na and Ring(d).norm((bx,by))==nb,'printed norms')
require(len(pairs)==40,'40 parsed pairs')
coeffsec=TEX.split('Coefficient bounds for the forty fixed pairs.')[1].split('\\end{longtable}')[0]
fixed={}
for m in re.finditer(r'^(\d+)&(\d+)&(\d+)&(\d+)&(\d+)&(\d+)',coeffsec,re.M):
    i,*v=map(int,m.groups());fixed[i]=v

def inventory():
    counts=[];triples=[]
    for d in [1,2,3,5,6,7,11,15,19,23]:
        R=Ring(d);found=set()
        for a in R.ball(2):
            na=R.norm(a)
            if na not in [1,2]:continue
            for r in R.ball(8):
                bn=R.mul(R.conj(a),R.sub(R.mul(r,r),(1,0)))
                if any(v%na for v in bn):continue
                b=tuple(v//na for v in bn);nb=R.norm(b)
                if not(9<=nb<=(63 if na==1 else 15)):continue
                if (a,b) in [((1,0),(3,0)),((-1,0),(-3,0))]:continue
                found.add((a,b))
        canonical={R.canonical(a,b) for a,b in found}
        expected={(a,b) for _,(dd,a,b) in pairs.items() if dd==d}
        require(canonical==expected,f'pair inventory d={d} diff={canonical^expected}')
        counts.append((d,len(canonical),len(found)))
    for i,(d,a,b) in pairs.items():
        R=Ring(d)
        values=[]
        squares={R.mul(z,z) for z in R.ball(120)}
        for s in R.ball(224):
            if R.norm(s)<R.norm(b) or s in (a,b):continue
            aa=R.add(R.mul(a,s),(1,0));bb=R.add(R.mul(b,s),(1,0))
            direct=aa in squares and bb in squares
            require(direct==(R.square(aa) and R.square(bb)),f'square test disagreement {i,s}')
            if direct:values.append(s)
        byroot=set()
        na=R.norm(a)
        for rt in R.ball(23):
            sn=R.mul(R.conj(a),R.sub(R.mul(rt,rt),(1,0)))
            if any(z%na for z in sn):continue
            s=tuple(z//na for z in sn)
            if R.norm(b)<=R.norm(s)<=224 and s not in (a,b) and R.square(R.add(R.mul(b,s),(1,0))):byroot.add(s)
        require(set(values)==byroot,f'two independent triple enumerations {i}')
        triples.extend((i,*s) for s in values)
    expected=[]
    sec=TEX.split('All remaining triples and the common cutoff procedure.')[1].split('\\end{longtable}')[0]
    for m in re.finditer(r'^(\d+) & (-?\d+) & (-?\d+) & (\d+) & (\d+) & (\d+) \\\\',sec,re.M):
        i,x,y,ns,mb,n0=m.groups();expected.append((int(i),int(x),int(y)))
    require(set(triples)-set(expected)=={(7,3,0)},f'triple difference {set(triples)^set(expected)}')
    require(len(triples)==57 and len(expected)==56,'triple counts')
    print('PAIR COUNTS',counts)
    print('PAIR TOTAL',sum(v for _,_,v in counts),'TRIPLES',len(triples),'MATCH 56+{(1,-3,3)}')
    return triples,counts

families=[('N2',2,2,2,F(5,4)),('N3',3,3,4,F(35,32)),('A3',3,2,4,F(35,16)),('A4',4,3,2,F(63,8)),('A5',5,4,4,F(495,64))]
def binom(a,n):
    v=F(1)
    for j in range(n): v=v*(a-j)/(j+1)
    return v
def make_polys():
    result=[]
    for name,k,m,g,Rold in families:
        pol=[];p=[m,k,k]
        for i in range(3):
            row=[]
            for j in [1,2]:
                out={};tot=p[j]+int(i==j)-1
                for h0 in range(tot+1):
                    for h1 in range(tot-h0+1):
                        h=[h0,h1,tot-h0-h1]
                        bij=4**k*binom(F(2*k+1,2),h[j])
                        for l in range(3):
                            if l!=j:bij*=binom(-p[l]-int(i==l),h[l])
                        if j==1:
                            bx=m;by=int(i==2)+h[2];pc=k-h[1];pd=k-int(i==2)-h[2]
                            bij*=(-1)**(k+int(i==2)+h[2])
                        else:
                            bx=int(i==1)+h[1];by=m;pc=k-h[2];pd=k-int(i==1)-h[1]
                        if i==0:pd-=1
                        if i==1:bx-=1
                        if i==2:by-=1
                        if j==1:bx-=1
                        if j==2:by-=1
                        require(min(bx,by,pc,pd)>=0,'termwise integral division')
                        for l in range(pc+1):
                            for w in range(pd+1):
                                ex=(bx+w+(l if j==1 else 0),by+pd-w+(l if j==2 else 0),l)
                                co=bij*comb(pc,l)*comb(pd,w)*(-1)**(pd-w)
                                out[ex]=out.get(ex,F(0))+co
                require(all(v.denominator==1 for v in out.values()),f'integral {name} {i,j}')
                out={ex:int(v) for ex,v in out.items() if v}
                require(all(r>=0 and s>=0 and 0<=l<=k and r+s==l+m+k-2 for r,s,l in out),'monomial weights')
                row.append(out)
            pol.append(row)
        print('POLYNOMIAL GENERATED',name,flush=True)
        result.append(pol)
    return result

POLYS=make_polys()

def evalpoly_byc(R,poly,a,b,k):
    vals=[(0,0) for _ in range(k+1)]
    for (r,s,l),n in poly.items():
        z=R.mul(R.power(a,r),R.power(b,s));vals[l]=R.add(vals[l],(n*z[0],n*z[1]))
    return vals
def fixed_coeff(R,a,b):
    nb=R.norm(b);output=[]
    for (name,k,m,g,Rold),rows in zip(families,POLYS):
        P=m+2*k;rs=[]
        for row in rows:
            ss=0
            for poly in row:
                for l,z in enumerate(evalpoly_byc(R,poly,a,b,k)):
                    r=k-l
                    if nb<=25:ss+=ceilroot(F(R.norm(z),nb**(P-2)),10)/1000**r
                    else:ss+=ceilroot(F(R.norm(z),nb**(P-2+r)),10)/200**r
            rs.append(F(101,100*g)*ss)
        output.append(int(max(rs))+1)
    return output

eps=F(1,4000);h=F(4001,3999);mu=F(51051,400000);kap=F(3535,409600)
Mm=F(8,9)*608**2;Rm=F(3920,81)
def adaptive(x):
    L=[5*x+10*x*x+10*x**4+5*x**5,10*(x+x**4+x**5),15+80*x+150*x*x+120*x**3+35*x**4,15+80*x+150*x*x+280*x**3+35*x**4]
    W=[295*x+143*x*x+143*x**4+295*x**5,590*(x+x**4+x**5),885+1144*x+435*x*x+42*x**3+170*x**4,885+1144*x+435*x*x+2634*x**3+170*x**4]
    return [l/256+F(21,1024*999)*w for l,w in zip(L,W)]
def quotients(n,b2,s2,G,rho,M,R,Rp,start=1024,second=256):
    M2,Mn,M3,M4,M5=M;R2,Rn,R3,R4,R5=R
    return [
      F(start,12100)*n**5*b2**5*(1-eps),
      F(second,36)*n**4*s2**2/(M2*M2*h),
      G**2*n**10*b2**2/(16*R2**2*2400**2*h),
      G**2*n**10*s2/(16*R2**2*2400**2*h),
      G**2*n**10/(16*Rp**2*Mn**2*h),
      G**2*n**10/(16*Rp**2*M3**2*h),
      G**4*n**18/(16*R3**2*Mm**2*h),
      G**4*n**17*b2*s2/(16*Rm**2*Mn**2*h),
      G**2*n**11*b2**2/(16*Rn**2*M4**2*h),
      G**4*n**18*b2**4/(16*Rn**2*4143**2*h),
      rho**2*n**7/(16*M4**2*h),
      G**4*n**19*b2/(144*R4**2*42**8*h),
      G**4*n**18/(4*M5**2*F(5,2)**8*h),
      G**6*n**27*b2**5/(144*R5**2*42**10*h),
      G**6*n**26*b2**2/(36*42**12*F(5,2)**10*h)]

def regions():
    regionsec=TEX.split('\\label{tab:pr-five-regions}')[1].split('\\end{tabular}')[0]
    coeffsec=TEX.split('\\label{tab:pr-five-coefficients}')[1].split('\\end{tabular}')[0]
    regionrows=[];coefrows=[]
    for line in regionsec.splitlines():
        if re.match(r'^(?:[1-4]|5--8)&',line):
            vals=line.rstrip('\\').split('&');regionrows.append([F(v) for v in vals[1:]])
    for line in coeffsec.splitlines():
        if re.match(r'^(?:[1-4]|5--8)&',line):
            vals=line.rstrip('\\').split('&');coefrows.append(([int(v) for v in vals[1:6]],F(vals[6])))
    require(len(regionrows)==len(coefrows)==5,'five current source regions')
    data=[(*row,*co) for row,co in zip(regionrows,coefrows)]
    output=[]
    for m,n,U,t,q,Fp,G,rho,Rp,M,margin in data:
        require(q==m/U and t*n*m*m>=54,'region ratio and retained-gap guard')
        l=F(0);r=F(2)
        for _ in range(30):
            z=(l+r)/2;val=5*z*z/q+4*(1+q**-2)*z-4/q
            if val<0:l=z
            else:r=z
        Fbound=(4-l*l)**4*(1+q**-2+r/q)**2
        require(Fbound<Fp,'phase bound')
        Q0=2*(1+1/q)**2;H=F(9999,10000)*F(200,19)**2*t**5*(1-mu*Q0/m**2)**2
        require(H/Q0**2>G and H/(kap*Fp)>rho,'G rho bounds')
        Mcalc=[]
        for (name,k,mm,g,Rold),rows in zip(families,POLYS):
            bound=F(101,100*g)*max(sum(F(abs(v))*q**(-r)*F(1000)**(l-k) for poly in row for (r,s,l),v in poly.items()) for row in rows)
            Mcalc.append(int(bound)+1)
        require(Mcalc==M,f'region M {n}: {Mcalc}!={M}')
        require(max(adaptive(1/q))<Rp,'adaptive region bound')
        R=[old*F(101,200)*((1+1/q)/2)**(2*k) for name,k,mm,g,old in families]
        qs=quotients(n,m*m,m*m,F(G),F(rho),M,R,Rp)
        require(min(qs)>margin,f'region margin {n}: {float(min(qs))}')
        require(F(G,105)*n**4*m**4>1,'tail')
        output.append({'norm_lower':int(n),'M':Mcalc,'quotients':[str(v) for v in qs],'minimum':str(min(qs)),'phase_upper':str(Fbound),'G_actual':str(H/Q0**2),'rho_actual':str(H/(kap*Fp))})
        print('REGION',n,'M',Mcalc,'min quotient',float(min(qs)),'junction',qs.index(min(qs))+1,'Fp',float(Fbound),'Gactual',float(H/Q0**2),'rhoactual',float(H/(kap*Fp)))

    return output

def fourth_gap():
    edges=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    radii=dict(zip(edges,[F(1,3),F(1,15),F(1,90),F(1,45),F(1,270),F(1,1350)]))
    ks=dict(zip(edges,[F(x,10000) for x in (5506,5087,5014,5029,5005,5001)]))
    k=lambda i,j:ks[tuple(sorted((i,j)))]
    for e in edges:require(F(1,2)<ks[e]<1 and radii[e]*ks[e]**2<2*ks[e]-1,'six rationalization inequalities')
    # Expand det(J+E) in six variables, directly by the 24 permutations.
    poly={}
    for p in permutations(range(4)):
        sign=(-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
        terms={(0,)*6:sign}
        for i,j in enumerate(p):
            if i==j:continue
            idx=edges.index(tuple(sorted((i,j))));nextterms=dict(terms)
            for ex,v in terms.items():
                ne=list(ex);ne[idx]+=1;ne=tuple(ne);nextterms[ne]=nextterms.get(ne,0)+v
            terms=nextterms
        for ex,v in terms.items():poly[ex]=poly.get(ex,0)+v
    poly={ex:v for ex,v in poly.items() if v}
    require(all(sum(ex) in (3,4) for ex in poly),'determinant has only cubic/quartic terms')
    require(sum(1 for ex in poly if sum(ex)==3)==22,'four triangles+twelve paths+six repeated-edge cubic monomials')
    # Leading cubic normalization: multiply by ab s t after setting Eij=1/(2ai aj).
    lead={}
    for ex,v in poly.items():
        if sum(ex)!=3:continue
        power=[1]*4
        for (i,j),e in zip(edges,ex):power[i]-=e;power[j]-=e
        key=tuple(power);lead[key]=lead.get(key,F(0))+F(v,8)
    expected={tuple(1 if i==j else -1 for i in range(4)):F(1,4) for j in range(4)}
    expected.update({tuple(-1 if i in e else 0 for i in range(4)):F(-1,4) for e in edges})
    lead={ex:v for ex,v in lead.items() if v}
    require(lead==expected,'exact signed leading cubic')
    W=F(1,4050);alpha=2*k(0,1)*k(0,2)*k(1,2)
    beta=6*k(1,2)*k(1,3)*k(2,3)*W+2*k(0,2)*k(0,3)*k(2,3)/450+2*k(0,1)*k(0,3)*k(1,3)/18
    defects=[]
    for i,j in edges:
        ell,m=[z for z in range(4) if z not in (i,j)]
        cc=2*(k(i,j)*(k(i,ell)*k(j,m)+k(i,m)*k(j,ell))+k(i,j)**2*k(ell,m))-F(1,2)
        require(cc>0,'all defect coefficients positive')
        defects.append(str(cc));beta+=cc*radii[i,j]
    beta+=(k(0,1)*k(2,3)+k(0,2)*k(1,3)+k(0,3)*k(1,2))**2*W
    require(alpha<F(141,500) and beta<F(21,125),'common tau constants')
    require(F(141,500)*F(59,20)+F(21,125)==F(9999,10000),'strict tau contradiction')
    return dict(alpha=str(alpha),beta=str(beta),defects=defects,determinant_terms=[dict(exponents=ex,coefficient=v) for ex,v in poly.items()])

def fixed41():
    output=[];tau=F(59,20)
    for i,d,a,b in [(i,*row) for i,row in pairs.items()]+[('rational',1,(1,0),(3,0))]:
        R=Ring(d);na=R.norm(a);nb=R.norm(b);nab=R.norm(R.sub(a,b))
        M=fixed_coeff(R,a,b)
        printed=fixed[i] if i!='rational' else [12,41,25,300,1100]
        require(M==printed if i!='rational' else all(x<=y for x,y in zip(M,printed)),f'actual source coefficient tuple {i}')
        A1=ceilroot(na,8);B0=floorroot(nb,8);B1=ceilroot(nb,8)
        Q0=ceilroot(F(nab,nb),8)*(1+A1/15)*(1+B1/15)
        require(1/Q0>mu/nb,'positive fixed-pair gap coefficient')
        K=F(9999,100)*na*(1/Q0-mu/nb)**2;G=K*tau**5
        require(K>F(7,4),'K>7/4 each41')
        H=F(9999,100)*na*tau**5*(1-mu*Q0/nb)**2
        F0=F(nab,nb)*(1+A1/15)**2*(1+B1/15)**8;rho=H/(kap*F0)
        require(B0/A1>2,'adaptive pair ratio>2')
        Rs=[old*F(101,200)*F(nab,4*nb)**k for _,k,_,_,old in families]
        Rp=max(adaptive(A1/B0));require(Rp==adaptive(A1/B0)[3],'combined V dominates residual')
        qs=quotients(1,nb,225,G,rho,printed,Rs,Rp,start=9604,second=784)
        require(min(qs)>F(59,50),'every41 pair every15 margin')
        tail=G*B0*15**3/105;require(tail>1,'every41 infinite tail')
        output.append(dict(pair=i,coefficient_minima=M,printed=printed,K=str(K),G=str(G),rho=str(rho),quotients=[str(q) for q in qs],minimum=str(min(qs)),tail=str(tail)))
    return output

def least(c,p):
    low=0;high=1
    while c*high**p<=1:high*=2
    while high-low>1:
        mid=(low+high)//2
        if c*mid**p>1:high=mid
        else:low=mid
    require(c*high**p>1 and c*(high-1)**p<=1,'minimal strict threshold')
    return high

def cutoff(R,a,b,s,M):
    na,nb,ns=map(R.norm,(a,b,s));dab=R.norm(R.sub(a,b));das=R.norm(R.sub(a,s));dbs=R.norm(R.sub(b,s))
    q1=ceilroot(F(dab*das*dbs,nb*ns*ns),10)
    require(1/q1>mu/nb,'fixed-triple positive difference')
    K=F(9999,100)*na*(1/q1-mu/nb)**2;require(K>1,'fixed-triple K>1')
    phi=ceilroot(max(F(dbs**6,ns**6),F(nb*dbs**5,ns**6)),10);Rbs=kap*phi
    M2,Mn,M3,M4,M5=M;Mbs=fixed_coeff(R,b,s)[1]
    R2,Rn,R3,R4,R5=[old*F(101,200)*F(dab,4*nb)**k for _,k,_,_,old in families]
    Rp=max(adaptive(1/floorroot(F(nb,na),10)))
    cs=[((1-eps)*K**4/(12100*ns**5),5),
      (K**2/(36*M2**2*nb**4*ns**2*h),4),
      (K**2/(16*R2**2*2400**2*nb**3*ns**5*h),5),
      (K**2/(16*R2**2*2400**2*nb**5*ns**4*h),5),
      (K**2/(16*Rp**2*M3**2*nb**6*ns**4*h),5),
      (K**4/(16*R3**2*Mm**2*nb**6*ns**12*h),10),
      (K**4/(16*Rm**2*Mn**2*nb**9*ns**9*h),10),
      (K**4/(16*Rn**2*Mbs**2*nb**3*ns**11*h),10),
      (K**2/(16*Rbs**2*M4**2*nb**5*ns**5*h),5),
      (K**4/(144*R4**2*42**8*nb**15*ns**4*h),10),
      (F(2,5)**8*K**4/(4*M5**2*nb**16*ns**4*h),10),
      (K**6/(144*R5**2*42**10*nb**19*ns**6*h),15),
      (F(2,5)**10*K**6/(36*42**12*nb**22*ns**6*h),15)]
    thresholds=[least(c,p) for c,p in cs]
    ntail=least(K*K/(105**2*nb**4*ns**2),5)
    eta=1 if na==1 else 4*na*nb*ns+1
    N0=max(1024,ntail,*thresholds,eta)
    return dict(q1=str(q1),K=str(K),Mbs=Mbs,thresholds=thresholds,ntail=ntail,eta=eta,N0=N0,quotients=[str(c*N0**p) for c,p in cs])

def complement(R,b,s,ns,N0):
    nb=R.norm(b);tail=4*nb*ns
    if tail>=N0:return dict(method='tail',tail=tail)
    q=R.add(b,s);hh=R.add((1,0),R.mul((2,0),R.mul(b,s)));delta=R.sub(R.power(R.sub(b,s),2),(4,0))
    nq,nh,nd=map(R.norm,(q,hh,delta))
    A=2*(floorroot(nh,10)-ceilroot(nq,10));B=ceilroot(nd,10)+2*ceilroot(nq,10)+1
    require(A>0,'quadratic lower linear coefficient')
    bound=F(N0)-A*floorroot(N0,10)+B
    if bound<-100:return dict(method='quadratic',tail=tail,norms=[nq,nh,nd],upper=str(bound))
    require((R.d,b,ns,N0)==(2,(-3,0),12,1024),'exact exceptional pair only')
    require(s in [(-2,-2),(-2,2)],'exact exceptional third only')
    require((nq,nh,nd)==(33,457,153),'exceptional norms')
    require(F(213,10)**2<457 and F(29,5)**2>33 and F(62,5)**2>153,'exceptional exact radical bounds')
    require(901-60*F(213,10)+62*F(29,5)+F(62,5)==-5,'f(30)<-5')
    residues={(u*u+2*v*v)%16 for u in range(16) for v in range(16)}
    require(14 not in residues and 15 not in residues,'annulus30,31 mod16')
    require(all(u%4==0 and v%4==0 for u in range(16) for v in range(16) if (u*u+2*v*v)%16==0),'annulus32 divisibility by4')
    solutions=[(u,v) for u in range(-5,6) for v in range(-4,5) if u*u+2*v*v==32]
    require(solutions==[(0,-4),(0,4)],'annulus32 root list')
    return dict(method='quadratic+annulus',tail=tail,norms=[nq,nh,nd],f30_upper='-5',residues=sorted(residues),solutions32=solutions)


def compact(s):
    return re.sub(r'\s+', '', s).replace('&','')

def source_field_controls():
    """Check actual printed formulas, then their algebra in polynomial rings."""
    normalized=compact(TEX)
    expected=[
        r'\omega^2-\tau_K\omega+\nu_K=0',
        r'(\tau_K,\nu_K)=(1,(d+1)/4)',
        r'\Delta_K=4\nu_K-\tau_K^2',
        r'(x,y)(v,w)=(xv-\nu_Kyw,xw+yv+\tau_Kyw)',
        r'N(x,y)=x^2+\tau_Kxy+\nu_Ky^2',
        r'\overline{(x,y)}=(x+\tau_Ky,-y)',
        r'4N(x,y)=(2x+\tau_Ky)^2+\Delta_Ky^2',
        r'\sqrt{4L/\Delta_K}',r'\sqrt{4L-\Delta_Ky^2}',
        r'e\equiv\tau_Ky\pmod2',r'x=(e-\tau_Ky)/2',
        r'p=2z_x+\tau_Kz_y',r'q=z_y',
        r'e^2=2n+p',r'\Delta_Kf^2=2n-p',r'ef=q',
        r'e\equiv\tau_Kf\pmod2',r'((e-\tau_Kf)/2,f)',
    ]
    for text in expected:
        require(compact(text) in normalized,'printed field formula: '+text)

    class P:
        """Sparse Z[x,y,v,w,tau,nu] for exact generic identities."""
        def __init__(self,a):
            self.a={e:c for e,c in a.items() if c} if isinstance(a,dict) else ({(0,)*6:a} if a else {})
        def __add__(self,b):
            b=b if isinstance(b,P) else P(b); d=dict(self.a)
            for e,c in b.a.items():d[e]=d.get(e,0)+c
            return P(d)
        __radd__=__add__
        def __neg__(self):return P({e:-c for e,c in self.a.items()})
        def __sub__(self,b):return self+-aspoly(b)
        def __rsub__(self,b):return aspoly(b)+-self
        def __mul__(self,b):
            b=aspoly(b);d={}
            for e,c in self.a.items():
                for f,v in b.a.items():
                    k=tuple(a+b for a,b in zip(e,f));d[k]=d.get(k,0)+c*v
            return P(d)
        __rmul__=__mul__
        def __pow__(self,n):
            o=P(1)
            for _ in range(n):o=o*self
            return o
        def __eq__(self,b):return self.a==aspoly(b).a
    def aspoly(a):return a if isinstance(a,P) else P(a)
    x,y,v,w,tau,nu=[P({tuple(int(i==j) for i in range(6)):1}) for j in range(6)]
    def pmul(a,b):return (a[0]*b[0]-nu*a[1]*b[1],a[0]*b[1]+a[1]*b[0]+tau*a[1]*b[1])
    def pnorm(a):return a[0]**2+tau*a[0]*a[1]+nu*a[1]**2
    def pconj(a):return (a[0]+tau*a[1],-a[1])
    a=(x,y);b=(v,w);disc=4*nu-tau*tau
    require(pmul(a,pconj(a))==(pnorm(a),P(0)),'generic conjugation and norm')
    require(pconj(pconj(a))==a,'generic conjugation involution')
    require(pnorm(pmul(a,b))==pnorm(a)*pnorm(b),'generic multiplicative norm')
    require(4*pnorm(a)==(2*x+tau*y)**2+disc*y*y,'generic disk identity')
    z=pmul(a,a);n=pnorm(a);e=2*x+tau*y;f=y;p=2*z[0]+tau*z[1];q=z[1]
    require(e*e==2*n+p,'generic square e equation')
    require(disc*f*f==2*n-p,'generic square f equation')
    require(e*f==q,'generic square cross equation')
    require(pnorm(z)==n*n,'generic square norm equation')
    require(e-tau*f==2*x,'generic root reconstruction and parity')

    bounded=[]
    for d in range(1,66):
        if any(d%(k*k)==0 for k in range(2,isqrt(d)+1)):continue
        R=Ring(d)
        # A coordinate box rigorously covering N(z)<=225: |y|<=sqrt(900/D),
        # |x|<=(30+|y|)/2<=24. Use a larger fixed box to cross-check disk enumeration.
        box={(xx,yy) for xx in range(-32,33) for yy in range(-32,33) if R.norm((xx,yy))<=225}
        require(set(R.ball(225))==box,f'disk formula versus coordinate box d={d}')
        roots={R.mul(z,z) for z in R.ball(15)}
        for z in box:
            require(R.square(z)==(z in roots),f'square equations versus root enumeration d={d},z={z}')
        require(R.square((0,0)),'zero square included')
        bounded.append((d,len(box),len(roots)))
    return dict(printed_formulas_checked=expected,generic_polynomial_identities=9,
                bounded_square_tests=bounded,scope='Symbolic identities hold in Z[x,y,v,w,tau,nu]; finite cross-tests cover every squarefree d<=65 and every element of norm<=225.')

if __name__=='__main__':
    print('SOURCE',SOURCE_HASH,flush=True)
    report=dict(article_sha256=SOURCE_HASH, status='PASS', scope='Sections 8 and 9: exact enumeration, coefficient bounds, interval certificates and finite complement.')
    report['field_controls']=source_field_controls()
    triples,counts=inventory();report['triples_before_exclusion']=len(triples);report['pair_counts']=counts
    report['regions']=regions()
    report['fourth_gap']=fourth_gap()
    report['fixed41']=fixed41()
    table=TEX.split('All remaining triples and the common cutoff procedure.')[1].split('\\end{longtable}')[0]
    rows=[]
    for m in re.finditer(r'^(\d+) & (-?\d+) & (-?\d+) & (\d+) & (\d+) & (\d+) \\\\',table,re.M):
        i,x,y,ns,mbs,n0=map(int,m.groups());d,a,b=pairs[i];R=Ring(d);s=(x,y)
        require(R.norm(s)==ns,'source final third norm')
        data=cutoff(R,a,b,s,fixed[i]);require(data['Mbs']==mbs and data['N0']==n0,f'source printed Mbs/N0 {i,s}')
        data.update(pair=i,s=s,snorm=ns,complement=complement(R,b,s,ns,n0));rows.append(data)
    require(len(rows)==56,'six-column final table exactly56')
    report['fixed56']=rows
    methods={kind:sum(z['complement']['method']==kind for z in rows) for kind in ('tail','quadratic','quadratic+annulus')}
    require(list(methods.values())==[38,16,2],'new finite complement partition')
    report['completion_counts']=methods
    rational=cutoff(Ring(1),(1,0),(3,0),(8,0),[12,41,25,300,1100])
    require(rational['q1']=='1822916667/5000000000' and F(rational['K'])>744,'rational q1 K')
    require(rational['Mbs']==34 and rational['N0']==1024,'rational Mbs cutoff')
    require(rational['thresholds']==[3,23,96,101,108,179,131,9,34,76,123,29,438] and rational['ntail']==14,'all rational thresholds')
    require([s for s in range(4,15) if isqrt(s+1)**2==s+1 and isqrt(3*s+1)**2==3*s+1]==[8],'rational third8')
    require(32<2*3*8,'rational tail contradiction')
    report['rational']=rational
    report['minimum_41_overlap']=str(min(F(z['minimum']) for z in report['fixed41']))
    report['maximum_quadratic_upper']=str(max(F(z['complement']['upper']) for z in rows if z['complement']['method']=='quadratic'))
    require(hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_HASH,'source unchanged during run')
    report['checks']=COUNTER
    (OUT/'one_small_element.json').write_text(json.dumps(report,indent=2))
    print('PASS',COUNTER,'exact checks')
    print('41-pair minimum overlap',format(float(F(report['minimum_41_overlap'])),'.12f'),'> 59/50; exact value in JSON')
    print('All56 Mbs/N0 reproduced, maximumN0',max(z['N0'] for z in rows),'and rational N0',rational['N0'])
    print('Completion partition',methods,'maximum quadratic upper bound',format(float(F(report['maximum_quadratic_upper'])),'.12f'))
    print('Output: one_small_element.json')
