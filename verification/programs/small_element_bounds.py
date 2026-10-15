#!/usr/bin/env python3
"""Exact inequalities and polynomial identities for the small-element cases."""
import sys
if sys.flags.optimize:
    raise SystemExit("Run without -O or -OO: the verification requires assertions.")
from pathlib import Path
from fractions import Fraction as F
from math import prod
import hashlib,json,os,tempfile
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data/article.tex'
EXPECTED_SOURCE_SHA='9cb6ed2154056e158d2c420012e94ef7c397c69d66c0fd9301a1bb52197dc432'
assert hashlib.sha256(SRC.read_bytes()).hexdigest()==EXPECTED_SOURCE_SHA
if os.environ.get('D1_RESULTS_DIR'):
 OUT=Path(os.environ['D1_RESULTS_DIR']);OUT.mkdir(parents=True,exist_ok=True)
else:
 (ROOT/'reruns').mkdir(exist_ok=True)
 OUT=Path(tempfile.mkdtemp(prefix='enumeration-',dir=ROOT/'reruns'))
checks=[]
def check(label, condition, evidence):
 assert condition,(label,evidence)
 checks.append({'label':label,'evidence':str(evidence)})
eps=F(1,4000)
check('crossblock geometric denominator',1-F(1,9)-F(1,81)-F(1,729)>F(7,8),1-F(1,9)-F(1,81)-F(1,729))
check('crossblock final constant',F(49,324)>F(1,7),F(49,324))
check('order-one determinant factor',4**5*prod(F(1,2)-v for v in range(-1,4))==-1440,4**5*prod(F(1,2)-v for v in range(-1,4)))
coef1=(8192*(1+F(11,8000)))**2*F(25*1001,16000)
check('four-root coefficient 10300',coef1<10300**2,coef1/F(10300**2))
coef2=192**2*F(1000,999)**7
check('four-root residual 193',coef2<193**2,coef2/F(193**2))
check('six-V gap',32-1>F(121,4),F(31)-F(121,4))
check('six-V lower 101',100/(1-eps)<101,100/(1-eps))
check('six-V upper 2/5',4/(9*(1+eps))>F(2,5),4/(9*(1+eps)))
for a,b,label in [(1134**10,10**31,'small h'),((F(15,2))**10,10**9,'small r'),(2277**5,10**17,'small u'),(42**100,10**163,'unit h'),(F(5,2)**5,10**2,'unit r'),(253**100,10**241,'unit u')]:check(label,a<b,F(a,b))
g=lambda x:7*x+F(212,5)+F(4047,25)/(x-4)
check('two-small envelope coefficients',(F(4)*F(17,5)+8*F(31,10)-4*F(-9,10),2*F(17,5)*F(-9,10)+42*4)==(42,F(4047,25)),(42,F(4047,25)))
check('S large fourth guard',(F(3,8)*10**8)**2>32*3,F(3,8)**2*10**16)
check('S large DC guard',9*10**16>45248*3,F(9*10**16,45248*3))
check('S lower constant',F(2,5)*F(3,8)**5==F(243,81920)>F(1,400),F(243,81920))
check('S cutoff at 49',49-F(227,5)-F(4047,25)/45==F(1,375),F(1,375))
check('g8',g(8)<139,g(8))
check('g49',g(49)==389-F(1,375),g(49))
check('general-three upper 37',36*(1+eps)<37,36*(1+eps))
check('T nonzero contradiction at72',72-F(347,5)-F(4047,25)/68>0,72-F(347,5)-F(4047,25)/68)
check('T zero F3 bound',F(1,4)*F(1000,999)**3<F(1,3),F(1,4)*F(1000,999)**3)
check('T zero coefficient1250',F(2,5)*5**5==1250,1250)
check('g72',g(72)<549,g(72))
M=10**8; C=10**100
check('bounded ST first term',F(30*M*M,C)<F(1,2),F(30*M*M,C))
check('bounded ST first lower weakened',F(C**5,36*M*M)-F(1,C)>F(C**5,100*M*M),F(1,36)-F(1,100))
check('bounded ST direct-four guard',F(C**5,100*M*M)>max(1000*C,17*9*M*M*C**3),F(C*C,1700*9*M**4))
check('bounded ST second term',F(78*M**5,C)<F(1,2),F(78*M**5,C))
check('bounded ST second lower weakened',F(625*C**7,2401*M**6)-F(1,C)>F(C**7,100*M**6),F(625,2401)-F(1,100))
check('bounded ST ordinary lower factor',F(512,81)*608**4*3**6<10**15,F(512,81)*608**4*3**6)
check('bounded ST ordinary upper denominator',2*28**4*F(10,81)**2*3**4<10**7,2*28**4*F(10,81)**2*3**4)
check('bounded ST overlap',10**143*C**5< F(C**7,10**50),F(C*C,10**193))
check('bounded ST upper envelope',g(100)<700+60,g(100)-700)
check('bounded ST upper over envelope',F(C**11,10**87)>10**60*C**7,F(C**4,10**147))
check('descent lambda contradiction',28>F(676,25),F(28)-F(676,25))
check('descent additive26',(F(20,14))**2>2,F(100,49))
check('descent strict decrease33',F(33,5)+26<33,F(33,5)+26)
check('unit cubic hbound',F(18,80)==F(9,40)<F(1,4),F(9,40))
check('unit cubic constant',F(16,15)**2*F(81,80)**3<2,F(16,15)**2*F(81,80)**3)
check('unit opposite pair contradiction',3*3**4-4*3**2-4>0,3*3**4-4*3**2-4)
check('unit cubic iterated denominator',2**13<10**4,2**13)
check('unit continuous parameters',(4*F(241,100)+8*F(163,100)-4*F(-2,5)+F(2,5))==F(617,25),F(617,25))
check('unit continuous residual',2*F(241,100)*F(-2,5)+(F(617,25)-F(2,5))*F(203,100)==F(118401,2500),F(118401,2500))
check('unit cutoff s4',80-F(717,25)-F(118401,4925)>0,80-F(717,25)-F(118401,4925))
check('small anchor F3 bound',F(1000,999)**3<4,F(1000,999)**3)
check('small anchor rationalization',F(1,12*1000**2)>F(1,1000**3),F(1000,12))
check('small anchor contradiction E<=4V',F(270,1000)+F(27,1000**2)<1,F(270,1000)+F(27,1000**2))
check('small anchor nonzero E lower',F(4*1000**5,81)-F(1,1000)>F(1000**5,100),F(4,81)-F(1,100))
check('small anchor final growth',F(1000**5,100)>5*1000**2,1000**3>500)
check('three-small envelope h',42**2*8**3<10**6,42**2*8**3)
check('three-small envelope u',253*8<10**4,253*8)
check('three-small envelope r',F(2,15)>F(1,100),F(2,15))
check('three-small envelope at40',48+F(224,35)+F(2,5)<60,48+F(224,35)+F(2,5))
check('three-small window lower',F(512,81)*608**4*3**22<10**40,F(512,81)*608**4*3**22)
check('three-small window upper',2*28**4*F(10,81)**2*3**14<10**22,2*28**4*F(10,81)**2*3**14)
check('four-fixed window lower',2*10300**2*3**12*85**12<10**38,2*10300**2*3**12*85**12)
check('four-fixed window upper',2*193**2*3**8*85**10<10**29,2*193**2*3**8*85**10)
check('three-large first lower',F(5,100**2)==F(1,2000),F(1,2000))
check('three-large second lower',F(5**5,100)>F(1,2000),F(5**5,100))
check('three-large lower versus upper',10**120>2000*10**60,10**120//(2000*10**60))
check('two-vanishing U bound',10**40>100*405**3,100*405**3)
check('four-anchor nonvanishing guard',10**80>32*3*3*85,32*3*3*85)
check('four-anchor coefficient',5*3*85**2==108375 and 3*85==255,(108375,255))
check('four-anchor first term',F(216750,10**40)<F(1,2),F(216750,10**40))
check('four-anchor lower W',F(10**200,260100)-F(1,10**40)>F(10**200,10**6),F(1,260100)-F(1,10**6))
check('four-anchor lower above window',10**80>10**44,10**80//10**44)
check('four-anchor windows overlap',10**80>10**69,10**80//10**69)
check('three-small upper contradiction',10**160>10**82,10**160//10**82)
# Sparse polynomial identities in the quotient by three monic quadratic relations.
# a,b,z,h,r,s with h^2=ab+1,r^2=az+1,s^2=bz+1.
def add(*ps):
 out={}
 for p in ps:
  for k,c in p.items():out[k]=out.get(k,0)+c
 return {k:c for k,c in out.items() if c}
def scale(p,c):return {k:v*c for k,v in p.items() if v*c}
def mul(p,q):
 out={}
 for k,c in p.items():
  for l,d in q.items():
   key=tuple(a+b for a,b in zip(k,l));out[key]=out.get(key,0)+c*d
 return out
def var(i):return {tuple(int(j==i) for j in range(6)):1}
one={(0,)*6:1};a,b,z,h,r,s=[var(i) for i in range(6)]
relations={3:add(mul(a,b),one),4:add(mul(a,z),one),5:add(mul(b,z),one)}
def reduce(p):
 out={};todo=list(p.items())
 while todo:
  k,c=todo.pop();i=next((i for i in (3,4,5) if k[i]>=2),None)
  if i is None:out[k]=out.get(k,0)+c;continue
  kk=list(k);kk[i]-=2
  for l,d in relations[i].items():todo.append((tuple(x+y for x,y in zip(kk,l)),c*d))
 return {k:c for k,c in out.items() if c}
zplus=add(a,b,mul(add(scale(mul(a,b),2),one),z),scale(mul(h,mul(r,s)),2))
rplus=add(mul(h,r),mul(a,s));splus=add(mul(b,r),mul(h,s))
for name,p in [('mutation az+1',add(mul(rplus,rplus),scale(add(mul(a,zplus),one),-1))),('mutation bz+1',add(mul(splus,splus),scale(add(mul(b,zplus),one),-1))),('mutation inverse',add(a,b,mul(add(scale(mul(a,b),2),one),zplus),scale(mul(h,mul(rplus,splus)),-2),scale(z,-1)))]:
 residual=reduce(p);check(name,not residual,residual)
# The opposite-pair factorization substitutes J^2=4a^2(1-b^2)(1-a^2b^2).
G=add(scale(mul(mul(a,a),mul(b,b)),2),scale(mul(a,a),-1),scale(one,-1))
J2=scale(mul(mul(a,a),mul(add(one,scale(mul(b,b),-1)),add(one,scale(mul(mul(a,a),mul(b,b)),-1)))),4)
right=mul(add(mul(a,a),scale(one,-1)),add(mul(a,a),scale(one,-1)))
res=add(mul(G,G),scale(J2,-1),scale(right,-1));check('opposite-pair factorization',not res,res)
assert hashlib.sha256(SRC.read_bytes()).hexdigest()==EXPECTED_SOURCE_SHA
report={'source_sha256':EXPECTED_SOURCE_SHA,'checks':checks,'outcome':'PASS','assertions':len(checks)}
(OUT/'small_element_bounds.json').write_text(json.dumps(report,indent=2))
print('Source SHA256',report['source_sha256'])
for row in checks:print('PASS',row['label'],row['evidence'])
print('PASS',len(checks),'exact checks')
