#!/usr/bin/env python3
"""Rebound v2 rerun of the independent v1 Section 7 checker. Computational core unchanged.
Doubled coordinates (u,v) mean (u+v sqrt(-d))/2. All decisions exact.
"""
from pathlib import Path
from math import isqrt
from fractions import Fraction as F
from itertools import combinations
import json,re,sys,time,hashlib
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'output/D1_bound5-v2.txt'
EXPECTED_SOURCE_SHA='5edcf8777371c4ad8979f122022d3119a8ae2de75d9674b336df994bba860b70'
assert hashlib.sha256(SRC.read_bytes()).hexdigest()==EXPECTED_SOURCE_SHA
TEXT=SRC.read_text()
OUT=ROOT/'revision_v2/checks'
ZERO=(0,0); ONE=(2,0)
def squarefree(d): return all(d%(p*p) for p in range(2,isqrt(d)+1))
def issq(n): return n>=0 and isqrt(n)**2==n
class Quadratic:
 def __init__(self,d): self.d=d
 def valid(self,z):
  u,v=z
  return (u-v)%2==0 if self.d%4==3 else u%2==v%2==0
 def add(self,x,y): return x[0]+y[0],x[1]+y[1]
 def neg(self,x): return -x[0],-x[1]
 def conj(self,x): return x[0],-x[1]
 def scale(self,x,k): return k*x[0],k*x[1]
 def mul(self,x,y):
  u,v=x; w,t=y
  a,b=u*w-self.d*v*t,u*t+v*w
  assert a%2==b%2==0
  return a//2,b//2
 def norm(self,x):
  u,v=x; n=u*u+self.d*v*v
  assert n%4==0
  return n//4
 def root(self,x):
  n=self.norm(x); h=isqrt(n)
  if h*h!=n:return None
  u,v=x; rr=2*h+u; ss=2*h-u
  if rr<0 or ss<0 or ss%self.d:return None
  r=isqrt(rr); s=isqrt(ss//self.d)
  if r*r!=rr or self.d*s*s!=ss or r*s!=abs(v):return None
  ans=(r,s if v>=0 else -s)
  if not self.valid(ans):return None
  assert self.mul(ans,ans)==x
  return ans
 def compatible(self,x,y):return self.root(self.add(self.mul(x,y),ONE)) is not None
 def ball(self,n):
  for v in range(-isqrt(4*n//self.d),isqrt(4*n//self.d)+1):
   a=isqrt(4*n-self.d*v*v)
   for u in range(-a,a+1):
    if self.valid((u,v)):yield(u,v)
 def canonical(self,vertices):
  return min(tuple(sorted(self.scale(self.conj(x) if c else x,s) for x in vertices)) for c in (0,1) for s in (-1,1))
 def initial(self,a,b):
  return {z for z in self.ball(1089) if self.compatible(a,z) and self.compatible(b,z)}
 def generate(self,a,b,exponent):
  p=self.mul(a,b); r0=self.root(self.add(p,ONE)); assert r0 not in (None,ZERO)
  cap=10**(2*exponent); values=self.initial(a,b); queue=sorted(values); j=0
  while j<len(queue):
   z=queue[j];j+=1
   r=self.root(self.add(self.mul(a,z),ONE));s=self.root(self.add(self.mul(b,z),ONE));assert r is not None and s is not None
   middle=self.add(self.add(a,b),self.mul(self.add(self.scale(p,2),ONE),z))
   last=self.scale(self.mul(r0,self.mul(r,s)),2)
   for w in (self.add(middle,last),self.add(middle,self.neg(last))):
    assert self.valid(w)
    if self.norm(w)<=cap and w not in values:
     # Explicit forward root identities certify every generated value.
     assert self.compatible(a,w) and self.compatible(b,w)
     values.add(w);queue.append(w)
  return values
 def modular_masks(self,p):
  # Integral coordinates x+y*w mod p, w^2=t*w-nu. p is odd.
  t=1 if self.d%4==3 else 0;nu=(self.d+1)//4 if t else self.d
  vals=[(x,y) for x in range(p) for y in range(p)]
  def mul(x,y):return ((x[0]*y[0]-nu*x[1]*y[1])%p,(x[0]*y[1]+x[1]*y[0]+t*x[1]*y[1])%p)
  squares={mul(x,x) for x in vals}; masks=[]
  for x in vals:
   mask=0
   for j,y in enumerate(vals):
    z=mul(x,y)
    if ((z[0]+1)%p,z[1]) in squares:mask|=1<<j
   masks.append(mask)
  return masks
 def graph(self,vertices):
  vertices=sorted(vertices); n=len(vertices); neighbors=[set() for _ in vertices]; checked=0;passed=0
  primes=(3,5,7,11,13)
  mods=[]
  for p in primes:
   inv2=(p+1)//2
   def residue(z):
    u,v=z
    if self.d%4==3:x,y=(u-v)//2,v
    else:x,y=u//2,v//2
    return (x%p)*p+y%p
   rs=[residue(z) for z in vertices]; masks=self.modular_masks(p)
   mods.append((rs,[masks[j] for j in rs],[1<<j for j in rs]))
  for i in range(n):
   rowmasks=[item[1][i] for item in mods]
   for j in range(i):
    if not all(rowmasks[k]&mods[k][2][j] for k in range(len(primes))):continue
    passed+=1
    if self.compatible(vertices[i],vertices[j]): neighbors[i].add(j);neighbors[j].add(i)
  edges=sum(map(len,neighbors))//2
  triangles=sum(len(neighbors[i]&neighbors[j]) for i in range(n) for j in neighbors[i] if j<i)//3
  return n,edges,triangles,passed

def pairs():
 rows=[]; allclasses=[]
 for d in range(1,36):
  if not squarefree(d):continue
  K=Quadratic(d);small=[z for z in K.ball(8) if z!=ZERO]
  classes=sorted({K.canonical((a,b)) for a,b in combinations(small,2) if K.compatible(a,b)})
  allclasses += [(d,a,b) for a,b in classes]
  for a,b in classes:
   if K.add(K.mul(a,b),ONE)==ZERO:continue
   started=time.monotonic(); values=K.generate(a,b,550)
   vertices={z for z in values if z not in (ZERO,a,b) and K.norm(z)>=9}
   counts=K.graph(vertices)
   row=[d,list(a),list(b),*counts];rows.append(row)
   print('PAIR',row,'seconds',round(time.monotonic()-started,2),flush=True)
 expected=[]
 section=TEXT.split('\\label{sm:bs7-pair-table}')[1].split('\\end{longtable}')[0]
 for m in re.finditer(r'(?m)^(\d+) & \$\((-?\d+),(-?\d+)\)\$ & \$\((-?\d+),(-?\d+)\)\$ & (\d+) & (\d+) ',section):
  d,u,v,x,y,n,e=map(int,m.groups());expected.append([d,[u,v],[x,y],n,e])
 actual=[r[:5] for r in rows]
 assert actual==expected,(actual,expected)
 assert len(allclasses)==72 and len(rows)==48
 assert sum(r[3] for r in rows)==36931 and sum(r[4] for r in rows)==36854 and not any(r[5] for r in rows)
 report={'hash':hashlib.sha256(SRC.read_bytes()).hexdigest(),'all_classes':allclasses,'rows':rows}
 (OUT/'small_v2_pairs.json').write_text(json.dumps(report,indent=2))
 print('PASS 72 classes, 48 nondegenerate graphs; 36931 vertices,36854 edges,zero triangles',flush=True)

def triples():
 listed=[(1,[(-4,-4),(-4,4),(4,0)],(0,2),(120,0),2),(1,[(-4,-2),(0,-2),(4,-2)],(0,1),(0,30),2),(1,[(-4,-2),(0,4),(4,-2)],(0,1),(0,-80),2),(1,[(-4,0),(-2,-4),(2,0)],(0,2),(8,24),4),(2,[(-4,0),(2,-2),(2,2)],(1,2),(-48,0),2),(3,[(-3,-1),(-3,1),(2,0)],(0,2),(16,0),2),(3,[(-3,-1),(1,-3),(2,0)],(0,2),(-24,8),4),(3,[(-1,-1),(0,-2),(1,-1)],(0,1),None,2),(7,[(-5,-1),(-5,1),(2,0)],(0,2),(48,0),2),(7,[(-5,-1),(0,2),(5,-1)],(0,1),(0,-64),2)]
 expected_counts=[61,89,64,105,71,112,112,57,83,55]; found={};quadruples=[]
 for d in range(1,36):
  if not squarefree(d):continue
  K=Quadratic(d);small=[z for z in K.ball(8) if z!=ZERO]
  adj={z:{w for w in small if w!=z and K.compatible(z,w)} for z in small}
  for tri in combinations(small,3):
   if all(y in adj[x] for x,y in combinations(tri,2)):
    can=K.canonical(tri);found.setdefault((d,can),set()).add(tuple(sorted(tri)))
    fourth=set.intersection(*(adj[z] for z in tri));assert not fourth,(d,tri,fourth)
 assert len(found)==10 and sum(map(len,found.values()))==24
 rows=[]
 for j,(d,tri,which,ext,orbit) in enumerate(listed):
  K=Quadratic(d);assert len(found[(d,K.canonical(tri))])==orbit
  a,b=[tri[i] for i in which]; values=K.generate(a,b,40)
  # Cross-check the root-first initial construction against direct norm-ball scanning.
  bound=isqrt(K.norm(a)*1089)+1;root_initial=set()
  for q in K.ball(bound):
   num=K.mul(K.conj(a),K.add(K.mul(q,q),(-2,0)));na=K.norm(a)
   if any(x%na for x in num):continue
   z=tuple(x//na for x in num)
   if K.valid(z) and K.norm(z)<=1089 and K.compatible(b,z):root_initial.add(z)
  assert root_initial==K.initial(a,b)
  actual={z for z in values if z not in [ZERO]+tri and all(K.compatible(x,z) for x in tri)}
  assert len(values)==expected_counts[j] and actual==({ext} if ext else set()),(j,len(values),actual)
  assert all(9<=K.norm(z)<85**2 for z in actual)
  rows.append({'row':j+1,'d':d,'orbit':orbit,'pair_count':len(values),'extensions':sorted(actual),'initial_count':len(root_initial)})
  print('TRIPLE',rows[-1],flush=True)
 (OUT/'small_v2_triples.json').write_text(json.dumps({'hash':hashlib.sha256(SRC.read_bytes()).hexdigest(),'rows':rows},indent=2))
 print('PASS 24 triples,10 orbits,no 4-clique; all ten radius1e40 pair closures and extension lists',flush=True)

def pell():
 R=10**550; umax=isqrt(R+1); rows=[];totals=[0,0,0,0,0]
 for d in range(7,10000,8):
  if not squarefree(d):continue
  totals[0]+=1;a0=isqrt(d);m=0;q=1;a=a0;uprev2=0;uprev1=1;vprev2=1;vprev1=0;verts=[];iterations=0
  while True:
   u=a*uprev1+uprev2;v=a*vprev1+vprev2;iterations+=1
   if u>umax:break
   if v>0 and u*u-d*v*v==2:verts +=[-(u*u-1),u*u-1]
   uprev2,uprev1=uprev1,u;vprev2,vprev1=vprev1,v
   m=q*a-m; numerator=d-m*m; assert numerator%q==0;q=numerator//q;a=(a0+m)//q
  # Necessary modular square filters, then exact rational square tests.
  primes=(16,9,5,7,11,13);sq={p:{x*x%p for x in range(p)} for p in primes}
  adj=[set() for _ in verts]
  residues=[[x%p for x in verts] for p in primes]
  for i,x in enumerate(verts):
   for j in range(i):
    y=verts[j]
    if (x>0)==(y>0):
     if any((residues[k][i]*residues[k][j]+1)%p not in sq[p] for k,p in enumerate(primes)):continue
     yes=issq(x*y+1)
    else:
     n=-(x*y+1);yes=n>=0 and n%d==0 and issq(n//d)
    if yes:adj[i].add(j);adj[j].add(i)
  edges=sum(map(len,adj))//2;triangles=sum(len(adj[i]&adj[j]) for i in range(len(verts)) for j in adj[i] if j<i)//3
  if verts:totals[1]+=1
  totals[2]+=len(verts);totals[3]+=edges;totals[4]+=triangles
  rows.append([d,len(verts),edges,triangles,iterations])
  if len(rows)%100==0:print('PELL progress',len(rows),totals,flush=True)
 assert totals==[1016,409,19612,9806,0],totals
 (OUT/'small_v2_pell.json').write_text(json.dumps({'hash':hashlib.sha256(SRC.read_bytes()).hexdigest(),'totals':totals,'rows':rows},indent=2))
 print('PASS Pell graphs',totals,flush=True)
if __name__=='__main__':
 print('Source SHA256',hashlib.sha256(SRC.read_bytes()).hexdigest(),flush=True)
 for arg in sys.argv[1:] or ['triples','pairs','pell']:
  started=time.perf_counter();passed=False
  try:
   globals()[arg]();passed=True
  finally:
   elapsed=time.perf_counter()-started
   current=hashlib.sha256(SRC.read_bytes()).hexdigest()
   assert current==EXPECTED_SOURCE_SHA,'Source changed during run'
   metadata={'mode':arg,'passed':passed,'seconds':elapsed,'source_sha256':current,'algorithm_origin':'Rebound independent v1 checker, same computational core'}
   (OUT/('small_v2_runtime_'+arg+'.json')).write_text(json.dumps(metadata,indent=2)+'\n')
   print('RUNTIME',json.dumps(metadata),flush=True)
