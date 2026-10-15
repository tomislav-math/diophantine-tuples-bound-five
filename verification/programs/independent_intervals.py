#!/usr/bin/env python3
"""Second exact implementation of interval, coefficient and cutoff checks.

Reconstructs five finite polynomial families independently, then checks the
all-large overlap chain, five one-small-entry regions, all 41 fixed-pair
certificates, and all 56 fixed-triple cutoffs plus (1,3,8). No result cache or
other verification program is imported.
"""
import sys
if sys.flags.optimize:
    raise SystemExit("Run without -O or -OO: the verification requires assertions.")
from fractions import Fraction as F
from math import comb, isqrt
from collections import defaultdict
import re, hashlib, json, os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'data/article.tex'
source_bytes=source.read_bytes()
EXPECTED='9cb6ed2154056e158d2c420012e94ef7c397c69d66c0fd9301a1bb52197dc432'
assert hashlib.sha256(source_bytes).hexdigest()==EXPECTED
s=source_bytes.decode('utf-8')
import tempfile
if os.environ.get('D1_RESULTS_DIR'):
    RESULTS = Path(os.environ['D1_RESULTS_DIR'])
    RESULTS.mkdir(parents=True, exist_ok=True)
else:
    (ROOT/'reruns').mkdir(exist_ok=True)
    RESULTS = Path(tempfile.mkdtemp(prefix=Path(__file__).stem+'-', dir=ROOT/'reruns'))
    print('Results directory:', RESULTS)
report=dict(status='PASS',article_sha256=EXPECTED,uniform_coefficients=[],regions=[],fixed_pair_coefficients=[],fixed_pair_overlaps=[],fixed_triple_cutoffs=[])

def binom(x,k):
    v=F(1)
    for j in range(k): v*=F(x-j,j+1)
    return v

def polynomials(k,m):
    p=[m,k,k]; out=[]
    for i in range(3):
        root=[]
        di=[int(i==j) for j in range(3)]
        for j in [1,2]:
            terms=defaultdict(F)
            nh=p[j]+di[j]-1
            for h0 in range(nh+1):
                for h1 in range(nh-h0+1):
                    h=[h0,h1,nh-h0-h1]
                    coeff=4**k*binom(F(2*k+1,2),h[j])
                    for ell in range(3):
                        if ell!=j: coeff*=binom(-p[ell]-di[ell],h[ell])
                    if j==1:
                        coeff*=(-1)**(k+di[2]+h[2]); r=m; t=di[2]+h[2]; n=k-h[1]; e=k-di[2]-h[2]
                    else: r=di[1]+h[1]; t=m; n=k-h[2]; e=k-di[1]-h[1]
                    if i==0: e-=1
                    elif i==1: r-=1
                    else: t-=1
                    if j==1: r-=1
                    else: t-=1
                    assert min(r,t,n,e)>=0
                    for ell in range(n+1):
                        for z in range(e+1):
                            rr=r+e-z+(ell if j==1 else 0); tt=t+z+(ell if j==2 else 0)
                            terms[(rr,tt,ell)]+=coeff*comb(n,ell)*comb(e,z)*(-1)**z
            root.append({key:v for key,v in terms.items() if v})
        out.append(root)
    return out

families=[('N2',2,2,2,F(5,4)),('N3',3,3,4,F(35,32)),('A3',3,2,4,F(35,16)),('A4',4,3,2,F(63,8)),('A5',5,4,4,F(495,64))]
polys={name:polynomials(k,m) for name,k,m,g,R in families}

def Mq(name,q):
    _,k,m,g,_=next(f for f in families if f[0]==name)
    return F(101,100*g)*max(sum(abs(v)*q**(-r)*F(1000)**(l-k) for col in row for (r,t,l),v in col.items()) for row in polys[name])

Muni=[]
for name,k,m,g,R in families:
    rows=[[sum(abs(v) for col in row for (r,t,l),v in col.items() if l==ell) for ell in range(k+1)] for row in polys[name]]
    print('polynomial sums',name,rows)
    val=Mq(name,F(1)); Muni.append(val.numerator//val.denominator+1)
    print('uniform exact',name,str(val),'floor+1',Muni[-1])
    report['uniform_coefficients'].append(dict(family=name,coefficient_sums=rows,bound=str(val),integer_bound=Muni[-1]))
assert Muni==[243,4143,2265,165763,3113013]

h=F(4001,3999); eps=F(1,4000); Mm=F(8,9)*608**2; Rm=F(3920,81)

def overlaps(n,m,S2,G,rho,Ms,Rs,Rp,Mbs=4143,first=(1024,256)):
    M2,Mn,M3,M4,M5=Ms; R2,Rn,R3,R4,R5=Rs; B2=F(m)**2
    return [F(first[0],12100)*n**5*B2**5*(1-eps),F(first[1],36)*n**4*S2**2/(M2**2*h),G**2*n**10*B2**2/(16*R2**2*2400**2*h),G**2*n**10*S2/(16*R2**2*2400**2*h),G**2*n**10/(16*Rp**2*Mn**2*h),G**2*n**10/(16*Rp**2*M3**2*h),G**4*n**18/(16*R3**2*Mm**2*h),G**4*n**17*B2*S2/(16*Rm**2*Mn**2*h),G**2*n**11*B2**2/(16*Rn**2*M4**2*h),G**4*n**18*B2**4/(16*Rn**2*Mbs**2*h),rho**2*n**7/(16*M4**2*h),G**4*n**19*B2/(144*R4**2*42**8*h),G**4*n**18/(4*M5**2*F(5,2)**8*h),G**6*n**27*B2**5/(144*R5**2*42**10*h),G**6*n**26*B2**2/(36*42**12*F(5,2)**10*h)]

def rp(q):
    x=1/q
    return F(15+80*x+150*x*x+280*x**3+35*x**4,256)+F(21,1024*999)*(885+1144*x+435*x*x+2634*x**3+170*x**4)

ru=[f[-1] for f in families]
lg=overlaps(F(9),F(3),F(9),F(400),400/F(35,32),Muni,ru,F(3))
print('all-large min',lg.index(min(lg))+1,str(min(lg)),float(min(lg)))
assert min(lg)==F(7834457305497600,5386911336360281)
assert min(lg)>F(7,5) and min(v for i,v in enumerate(lg) if i!=10)>F(13,10)
assert rp(F(1))<3
report['all_large']=dict(quotients=[str(x) for x in lg],minimum=str(min(lg)),residual=str(rp(F(1))))
mu=F(51051,400000); kap=F(3535,409600)
regions=[(1,8,1,1,F(27,8),F(8),269,7495,20687,F(7,50),[36,282,145,4605,38390],F(1260,1000)),(2,4,2,F(3,2),F(3),F(8,3),368,1771,7973,F(7,20),[69,697,326,14771,171856],F(2329,1000)),(3,3,3,F(7,4),F(11,4),F(12,7),540,598,3226,F(18,25),[109,1306,663,34919,482131],F(1166,1000)),(4,3,4,F(2),F(29,10),F(3,2),637,624,3507,F(47,50),[130,1652,863,48019,701908],F(5459,1000)),(5,3,5,F(17,6),F(3),F(18,17),1121,375,2217,F(201,100),[220,3553,1944,134712,2400382],F(1321,1000))]
for label,m,n,U,t,q,Fp,G,rho,Rp,Ms,margin in regions:
    Q0=2*(1+1/q)**2; H=F(9999,10000)*F(200,19)**2*t**5*(1-mu*Q0/m**2)**2
    assert H/Q0**2>G and H/(kap*Fp)>rho and rp(q)<Rp
    lo=F(0); hi=F(2)
    for _ in range(30):
        mid=(lo+hi)/2
        if 5*mid**2/q+4*(1+q**-2)*mid-4/q>0: hi=mid
        else: lo=mid
    assert (4-lo**2)**4*(1+q**-2+hi/q)**2<Fp
    assert [Mq(name,q).numerator//Mq(name,q).denominator+1 for name,*_ in families]==Ms
    Rs=[R*F(101,200)*((1+1/q)/2)**(2*k) for name,k,mm,g,R in families]
    ov=overlaps(F(n),F(m),F(m*m),F(G),F(rho),Ms,Rs,Rp)
    assert min(ov)>margin
    assert F(G,105)*n**4*m**4>1
    print('region',label,'minimum',float(min(ov)), 'junction',ov.index(min(ov))+1)
    report['regions'].append(dict(region=label,coefficients=Ms,quotients=[str(x) for x in ov],minimum=str(min(ov)),required_margin=str(margin)))

# Exact ring evaluation of polynomials, using the field basis in the manuscript.
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mult(a,b,d):
    t=1 if d%4==3 else 0; n=F(d+1,4) if t else d
    return (a[0]*b[0]-n*a[1]*b[1],a[0]*b[1]+a[1]*b[0]+t*a[1]*b[1])
def power(a,n,d):
    z=(F(1),F(0))
    for _ in range(n): z=mult(z,a,d)
    return z
def norm(a,d):
    t=1 if d%4==3 else 0; n=F(d+1,4) if t else d
    return a[0]**2+t*a[0]*a[1]+n*a[1]**2

def sqrtbox(v,prec,up):
    scale=10**prec; n=isqrt((v*scale**2).numerator//(v*scale**2).denominator)
    if up and F(n*n,scale*scale)<v: n+=1
    return F(n,scale)

def Mfixed(name,a,b,d):
    _,k,m,g,_=next(f for f in families if f[0]==name); P=m+2*k; nb=norm(b,d)
    sums=[]
    for row in polys[name]:
        z=F(0)
        for col in row:
            for ell in range(k+1):
                value=(F(0),F(0))
                for (r,t,l),coeff in col.items():
                    if l==ell:
                        ab=mult(power(a,r,d),power(b,t,d),d); value=add(value,(coeff*ab[0],coeff*ab[1]))
                np=norm(value,d); rr=k-ell
                if nb<=25: z+=sqrtbox(np/nb**(P-2),10,True)/1000**rr
                else: z+=sqrtbox(np/nb**(P-2+rr),10,True)/200**rr
        sums.append(z)
    value=max(sums)*F(101,100*g)
    return value.numerator//value.denominator+1

pairpat=r'^(\d+)&(\d+)&\$\((-?\d+),(-?\d+)\)\$&\$\((-?\d+),(-?\d+)\)\$&(\d+)&(\d+)\\\\'
pairs=[]
for match in re.finditer(pairpat,s,re.M):
    i,d,ax,ay,bx,by,na,nb=map(int,match.groups()); pairs.append((i,d,(F(ax),F(ay)),(F(bx),F(by)),na,nb))
coefblock=s.split('\\caption{Coefficient bounds for the forty fixed pairs.}')[1].split('\\end{longtable}')[0]
cm={int(m[0]):list(map(int,m[1:])) for m in re.findall(r'^(\d+)&(\d+)&(\d+)&(\d+)&(\d+)&(\d+)\\\\',coefblock,re.M)}
assert len(pairs)==40 and len(cm)==40
for i,d,a,b,na,nb in pairs:
    assert norm(a,d)==na and norm(b,d)==nb
    calc=[Mfixed(name,a,b,d) for name,*_ in families]
    assert calc==cm[i],(i,calc,cm[i])
    report['fixed_pair_coefficients'].append(dict(pair=i,coefficients=calc))
print('fixed pair coefficient rows reconstructed',len(pairs))
rational=(-1,1,(F(1),F(0)),(F(3),F(0)),1,9)
assert [Mfixed(name,rational[2],rational[3],1) for name,*_ in families]==[12,41,25,300,1100]
# Fixed-pair S>=15 overlaps.
mins=[]
for i,d,a,b,na,nb in pairs+[rational]:
    na=F(na); nb=F(nb); nab=norm((a[0]-b[0],a[1]-b[1]),d)
    A1=sqrtbox(na,8,True); B0=sqrtbox(nb,8,False); B1=sqrtbox(nb,8,True)
    d1=sqrtbox(nab/nb,8,True); Q0=d1*(1+A1/15)*(1+B1/15); tau=F(59,20)
    K=F(9999,100)*na*(1/Q0-mu/nb)**2; G=K*tau**5
    H=F(9999,100)*na*tau**5*(1-mu*Q0/nb)**2
    F0=(nab/nb)*(1+A1/15)**2*(1+B1/15)**8; rho=H/(kap*F0)
    assert 1/Q0>mu/nb and K>F(7,4)
    Rs=[R*F(101,200)*(nab/(4*nb))**k for name,k,mm,g,R in families]
    Ms=cm[i] if i>=0 else [12,41,25,300,1100]
    # overlaps accepts m but uses only m^2, so evaluate with exact square substitution in a separate factor map.
    ov=overlaps(F(1),F(1),F(225),G,rho,Ms,Rs,rp(B0/A1),first=(9604,784))
    b_powers=[5,0,2,0,0,0,0,1,2,4,0,1,0,5,2]
    ov=[v*nb**p for v,p in zip(ov,b_powers)]
    assert min(ov)>F(59,50),(i,float(min(ov)))
    assert G*B0*15**3/105>1
    mins.append((min(ov),i,ov.index(min(ov))+1))
    report['fixed_pair_overlaps'].append(dict(pair=i,quotients=[str(x) for x in ov],minimum=str(min(ov)),tail=str(G*B0*15**3/105)))
print('S15 minimum',float(min(mins)[0]), 'pair/junction',min(mins)[1:])
print('All 41 fixed-pair certificates pass')
# Independent reconstruction of all 57 common cutoffs.
def threshold(c,p):
    hi=1
    while c*hi**p<=1: hi*=2
    lo=0
    while lo+1<hi:
        mid=(lo+hi)//2
        if c*mid**p>1: hi=mid
        else: lo=mid
    assert c*hi**p>1 and c*(hi-1)**p<=1
    return hi

def cutoffs(d,a,b,z,Ms):
    aa=norm(a,d); bb=norm(b,d); ss=norm(z,d)
    dab=norm((a[0]-b[0],a[1]-b[1]),d); das=norm((a[0]-z[0],a[1]-z[1]),d); dbs=norm((b[0]-z[0],b[1]-z[1]),d)
    q1=sqrtbox(dab*das*dbs/(bb*ss**2),10,True)
    K=F(9999,100)*aa*(1/q1-mu/bb)**2
    assert 1/q1>mu/bb and K>1
    phi=sqrtbox(max(dbs**6/ss**6,bb*dbs**5/ss**6),10,True); Rbs=kap*phi
    Mbs=Mfixed('N3',b,z,d)
    M2,Mn,M3,M4,M5=Ms
    R2,Rn,R3,R4,R5=[R*F(101,200)*(dab/(4*bb))**k for name,k,mm,g,R in families]
    q=sqrtbox(bb/aa,10,False); assert q>=1 and q*q<=bb/aa
    Rp=rp(q)
    cs=[((1-eps)*K**4/(12100*ss**5),5), (K*K/(36*M2*M2*bb**4*ss**2*h),4), (K*K/(16*R2**2*2400**2*bb**3*ss**5*h),5), (K*K/(16*R2**2*2400**2*bb**5*ss**4*h),5), (K*K/(16*Rp**2*M3**2*bb**6*ss**4*h),5), (K**4/(16*R3**2*Mm**2*bb**6*ss**12*h),10), (K**4/(16*Rm**2*Mn**2*bb**9*ss**9*h),10), (K**4/(16*Rn**2*Mbs**2*bb**3*ss**11*h),10), (K*K/(16*Rbs**2*M4**2*bb**5*ss**5*h),5), (K**4/(144*R4**2*42**8*bb**15*ss**4*h),10), (F(2,5)**8*K**4/(4*M5**2*bb**16*ss**4*h),10), (K**6/(144*R5**2*42**10*bb**19*ss**6*h),15), (F(2,5)**10*K**6/(36*42**12*bb**22*ss**6*h),15)]
    ns=[threshold(c,p) for c,p in cs]; ni=threshold(K*K/(105**2*bb**4*ss**2),5)
    N0=max(1024,ni,*ns,1 if aa==1 else 4*aa*bb*ss+1)
    return Mbs,N0,ns,ni,q1,K
tripleblock=s.split('\\label{tab:fixed56}')[1].split('\\end{longtable}')[0]
triples=[list(map(int,row)) for row in re.findall(r'^(\d+)\s*&\s*(-?\d+)\s*&\s*(-?\d+)\s*&\s*(\d+)\s*&\s*(\d+)\s*&\s*(\d+)\s*\\\\',tripleblock,re.M)]
assert len(triples)==56
for idx,u,v,ss,Mbs,N0 in triples:
    i,d,a,b,na,nb=pairs[idx]; z=(F(u),F(v)); assert norm(z,d)==ss
    actual=cutoffs(d,a,b,z,cm[i]); assert actual[:2]==(Mbs,N0),(idx,u,v,actual[:2],Mbs,N0)
    report['fixed_triple_cutoffs'].append(dict(pair=idx,third=[u,v],Mbs=actual[0],N0=actual[1],thresholds=actual[2],tail_threshold=actual[3],q1=str(actual[4]),K=str(actual[5])))
print('fixed 56 cutoffs reconstructed',len(triples))
rational_c=cutoffs(1,(F(1),F(0)),(F(3),F(0)),(F(8),F(0)),[12,41,25,300,1100])
assert rational_c[:4]==(34,1024,[3,23,96,101,108,179,131,9,34,76,123,29,438],14)
assert rational_c[4]==F(1822916667,5000000000) and rational_c[5]>744
print('rational cutoff',rational_c[:4])
report['rational_triple']=dict(Mbs=rational_c[0],N0=rational_c[1],thresholds=rational_c[2],tail_threshold=rational_c[3],q1=str(rational_c[4]),K=str(rational_c[5]))
assert hashlib.sha256(source.read_bytes()).hexdigest()==EXPECTED
out=RESULTS/'independent_intervals.json'
out.write_text(json.dumps(report,indent=2,default=str)+'\n')
print('PASS: five families, five regions, 41 pairs, 56 fixed triples and (1,3,8)')
print('Output:',out.name)
