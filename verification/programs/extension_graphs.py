#!/usr/bin/env python3
"""Complete finite extension graphs for the small-element cases.
Doubled coordinates (u,v) mean (u+v sqrt(-d))/2. All decisions are exact.
Run with triples, pairs, or pell to select the computation.
"""
from pathlib import Path
from math import isqrt
from fractions import Fraction as F
from itertools import combinations
import json,re,sys,time,hashlib,os,tempfile,ast
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data/article.txt'
EXPECTED_SOURCE_SHA='9cb6ed2154056e158d2c420012e94ef7c397c69d66c0fd9301a1bb52197dc432'
assert hashlib.sha256(SRC.read_bytes()).hexdigest()==EXPECTED_SOURCE_SHA
if os.environ.get('D1_RESULTS_DIR'):
 OUT=Path(os.environ['D1_RESULTS_DIR']);OUT.mkdir(parents=True,exist_ok=True)
else:
 (ROOT/'reruns').mkdir(exist_ok=True)
 OUT=Path(tempfile.mkdtemp(prefix='enumeration-',dir=ROOT/'reruns'))
TEXT=SRC.read_text()
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
 print('PAIR columns: field_d, anchor_a, anchor_b, vertices, edges, triangles, modular_filter_survivors',flush=True)
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
 report={'outcome':'PASS','source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'radius_exponent':550,'columns':['field_d','anchor_a','anchor_b','vertices','edges','triangles','modular_filter_survivors'],'totals':{'classes':len(allclasses),'graphs':len(rows),'vertices':sum(r[3] for r in rows),'edges':sum(r[4] for r in rows),'triangles':sum(r[5] for r in rows)},'all_classes':allclasses,'rows':rows}
 (OUT/'small_pairs.json').write_text(json.dumps(report,indent=2))
 print('PASS 72 classes, 48 nondegenerate graphs; 36931 vertices,36854 edges,zero triangles',flush=True)

def triples():
 # Parse affine expressions in t from the printed table using a restricted AST.
 def coordinate(expr):
  expr=re.sub(r'(\d)t',r'\1*t',expr)
  def read(n):
   if isinstance(n,ast.Constant) and isinstance(n.value,int):return F(n.value),F(0)
   if isinstance(n,ast.Name) and n.id=='t':return F(0),F(1)
   if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):
    a,b=read(n.operand);return -a,-b
   if isinstance(n,ast.BinOp):
    a,b=read(n.left);c,d=read(n.right)
    if isinstance(n.op,ast.Add):return a+c,b+d
    if isinstance(n.op,ast.Sub):return a-c,b-d
    if isinstance(n.op,ast.Mult):
     assert b*d==0;return a*c,a*d+b*c
    if isinstance(n.op,ast.Div):
     assert d==0 and c!=0;return a/c,b/c
   raise ValueError('Unsupported table expression: '+expr)
  a,b=read(ast.parse(expr,mode='eval').body)
  assert (2*a).denominator==(2*b).denominator==1
  return int(2*a),int(2*b)
 section=TEXT.split('Small triple&Extension&Orbit size')[1].split('\\bottomrule')[0]
 listed=[]
 pair_choices=[(0,2),(0,1),(0,1),(0,2),(1,2),(0,2),(0,2),(0,1),(0,2),(0,1)]
 for line in section.splitlines():
  if not re.match(r'^\d+&\d+&',line):continue
  ds,_,tr,ex,orb=line.rstrip('\\').split('&')
  tr=tr.replace('$','').replace('\\{','').replace('\\}','')
  tri=[coordinate(x) for x in tr.split(',')]
  ext=None if ex=='none' else coordinate(ex.strip('$'))
  listed.append((int(ds),tri,pair_choices[len(listed)],ext,int(orb)))
 assert len(listed)==10
 # These two displayed lists determine the selected pairs and closure counts.
 pair_block=TEXT.split('In the successive ten rows choose the pairs')[1].split('Each chosen pair')[0]
 parsed_pairs=[(int(i)-1,int(j)-1) for i,j in re.findall(r'\((\d),(\d)\)',pair_block)]
 assert parsed_pairs==pair_choices
 closure_block=TEXT.split('The numbers of distinct pair extensions obtained, including zero and small entries, are respectively')[1].split('Testing the third anchor')[0]
 expected_counts=[int(x) for x in re.search(r'61,89,64,105,71,112,112,57,83,55',closure_block).group().split(',')]
 found={};quadruples=[]
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
 (OUT/'small_triples.json').write_text(json.dumps({'outcome':'PASS','source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'radius_exponent':40,'total_triples':24,'orbits':10,'four_cliques':0,'rows':rows},indent=2))
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
 expected_match=re.search(r'(?m)^1016 &409 &19612 &9806 &0',TEXT)
 assert expected_match,'Pell graph totals missing from the article'
 expected=[int(x.strip()) for x in expected_match.group().split('&')]
 assert totals==expected,totals
 (OUT/'pell_graphs.json').write_text(json.dumps({'outcome':'PASS','source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'radius_exponent':550,'total_columns':['fields','nonempty_graphs','vertices','edges','triangles'],'totals':totals,'row_columns':['field_d','vertices','edges','triangles','continued_fraction_iterations'],'rows':rows},indent=2))
 print('PASS Pell graphs',totals,flush=True)
if __name__=='__main__':
 print('Source SHA256',hashlib.sha256(SRC.read_bytes()).hexdigest(),flush=True)
 for arg in sys.argv[1:] or ['triples','pairs','pell']:
  assert arg in ('triples','pairs','pell'),'Expected triples, pairs, or pell'
  started=time.perf_counter();passed=False
  try:
   globals()[arg]();passed=True
  finally:
   elapsed=time.perf_counter()-started
   current=hashlib.sha256(SRC.read_bytes()).hexdigest()
   assert current==EXPECTED_SOURCE_SHA,'Source changed during run'
   metadata={'mode':arg,'passed':passed,'seconds':elapsed,'source_sha256':current,'algorithm':'exact integer square tests; reversible descent or continued fractions'}
   (OUT/('extension_graphs_runtime_'+arg+'.json')).write_text(json.dumps(metadata,indent=2)+'\n')
   print('RUNTIME',json.dumps(metadata),flush=True)
