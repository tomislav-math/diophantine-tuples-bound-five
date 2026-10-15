#!/usr/bin/env python3
"""Reconstruct the forbidden intervals and their 15 overlap bounds exactly.

Checks Section 5.5 and Proposition 6.1. Laurent monomials are compared
symbolically, with exact rational coefficients. No earlier article is read.
"""
import sys
if sys.flags.optimize:
    raise SystemExit("Run without -O or -OO: the verification requires assertions.")
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import os

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'data/article.tex'
candidate = SOURCE.read_text(encoding='utf-8')
EXPECTED = '9cb6ed2154056e158d2c420012e94ef7c397c69d66c0fd9301a1bb52197dc432'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED
import tempfile
if os.environ.get('D1_RESULTS_DIR'):
    RESULTS = Path(os.environ['D1_RESULTS_DIR'])
    RESULTS.mkdir(parents=True, exist_ok=True)
else:
    (ROOT/'reruns').mkdir(exist_ok=True)
    RESULTS = Path(tempfile.mkdtemp(prefix=Path(__file__).stem+'-', dir=ROOT/'reruns'))
    print('Results directory:', RESULTS)

def table_at(marker):
    pos=candidate.index(marker)
    begin=candidate.index(r'\begin{tabular}',pos)
    end=candidate.index(r'\end{tabular}',begin)+len(r'\end{tabular}')
    return candidate[begin:end]

new_shared=table_at(r'\label{tab:pr-windows}')
new_overlaps=table_at(r'\label{tab:pr-overlaps}')
# A Laurent monomial in A,B,S,T,C, with exact rational coefficient.
letters='ABSTC'
class Mon:
    def __init__(self,c=1,e=None):
        self.c=F(c); self.e=tuple(e or (0,)*5)
    def __mul__(self,b):
        return Mon(self.c*b.c,[x+y for x,y in zip(self.e,b.e)])
    def __truediv__(self,b):
        return Mon(self.c/b.c,[x-y for x,y in zip(self.e,b.e)])
    def __pow__(self,n):
        assert n.denominator==1
        n=n.numerator
        return Mon(self.c**n,[x*n for x in self.e])
    def __add__(self,b):
        assert self.e==b.e
        return Mon(self.c+b.c,self.e)
    def __sub__(self,b):
        assert self.e==b.e
        return Mon(self.c-b.c,self.e)
    def __eq__(self,b):
        return self.c==b.c and self.e==b.e
    def value(self,vals):
        ans=self.c
        for a,n in zip(letters,self.e):ans*=F(vals[a])**n
        return ans

uniform=dict(M2=243,Mn=4143,M3=2265,M4=165763,M5=3113013,Mbs=4143,
             Mm=F(8,9)*608**2,R2=F(5,4),Rn=F(35,32),R3=F(35,16),
             R4=F(63,8),R5=F(495,64),Rbs=F(35,32),Rp=3,Rm=F(3920,81),
             h=F(4001,3999),eps=F(1,4000),G=400,rho=F(400*32,35))

class Parser:
    def __init__(self,s,vals,j=4):
        s=s.replace(r'\(', '').replace(r'\)','').replace('$','')
        s=s.replace(r'\cdot','*').replace(r'\varepsilon','eps').replace(r'\rho','rho')
        s=s.replace(r'\,','*')
        s=re.sub(r'([MR])_\{(bs)\}',r'\1\2',s)
        s=re.sub(r'([MR])_(bs|[2n345mp])',r'\1\2',s)
        s=s.replace('{','(').replace('}',')')
        self.tokens=re.findall(r'(?:[MR](?:bs|[2n345mp]))|eps|rho|\d+|[ABSTCGhj]|[()+*/^\-]',s)
        assert ''.join(self.tokens)==re.sub(r'\s','',s),s
        self.i=0; self.vals=dict(vals,j=Mon(j))
    def peek(self): return self.tokens[self.i] if self.i<len(self.tokens) else None
    def pop(self): t=self.peek();self.i+=1;return t
    def atom(self):
        t=self.pop()
        if t=='(':
            a=self.add();assert self.pop()==')';return a
        if t.isdigit():return Mon(int(t))
        return self.vals[t]
    def power(self):
        a=self.atom()
        if self.peek()=='^':
            self.pop();n=self.atom();assert not any(n.e)
            a=a**n.c
        return a
    def prod(self):
        a=self.power()
        while self.peek() not in (None,')','+','-'):
            t=self.peek()
            if t in ('*','/'):
                self.pop();b=self.power();a=a*b if t=='*' else a/b
            else:a=a*self.power()
        return a
    def add(self):
        a=self.prod()
        while self.peek() in ('+','-'):
            t=self.pop();b=self.prod();a=a+b if t=='+' else a-b
        return a
    def parse(self):
        a=self.add();assert self.peek() is None;return a

def formula(s,branch=0,j=4,unit_A=False):
    s=s.replace(r'\max(BS^3,B^3S^2)', ['(BS^3)','(B^3S^2)'][branch])
    s=s.replace('H_P', ['(BS^3)','(B^3S^2)'][branch])
    v={k:Mon(x) for k,x in uniform.items()}
    for k,a in enumerate(letters):
        e=[0]*5;e[k]=1;v[a]=Mon(1,e)
    if unit_A:v['A']=Mon(1)
    return Parser(s,v,j).parse()

def rows(tab):
    out={}
    for line in tab.splitlines():
        if '&' not in line:continue
        cells=line.rstrip().removesuffix(r'\\').split('&')
        label=cells[0].replace(r'\(','').replace(r'\)','').replace('$','').strip()
        if label not in ['Q','N_2','P','A_3','M','N_3','N_{3,bs}','A_4','A_5','O_j']:continue
        out[label]=[s.strip() for s in cells[1:3]]
    assert len(out)==10
    return out


# For 1 < U/Z + V Z, Z^2 = |cd+1|, the raw interval endpoints are
# l=4 U^2/C and u=1/(4 V^2 C). The factors 1/(1-epsilon),1/(1+epsilon)
# are applied separately as in the article. U,V below are the preceding
# approximation inequalities, not another copy of the endpoint table.
bound_terms={
    'Q': ('55T^3C^2','3T/(2A^2C^4)'),
    'N_2': ('M_2B^4C^3','R_2B^2/(A^3C^4)'),
    'P': ('2400H_PC^3','R_pS^2/(A^3C^5)'),
    'A_3': ('M_3B^6C^4','R_3B^3/(A^4C^5)'),
    'M': ('M_mB^3S^8C^3','R_mB^2S^5/(A^3C^6)'),
    'N_3': ('M_nB^7C^4','R_nB^3/(A^4C^6)'),
    'N_{3,bs}': ('M_{bs}S^7C^4','R_{bs}S^3/(B^4C^6)'),
    'A_4': ('M_4B^9C^5','R_4B^4/(A^5C^7)'),
    'A_5': ('M_5B^{12}C^6','R_5B^5/(A^6C^9)'),
    'O_j': ('3(42B^3C)^jC/B','(5B/(2AC^2))^j/2'),
}
common_rows=rows(new_shared)
comparisons=[]
for family,(u,v) in bound_terms.items():
    for branch in (0,1) if family=='P' else (0,):
        for j in (4,5,6) if family=='O_j' else (4,):
            U,V=formula(u,branch,j),formula(v,branch,j)
            reconstructed=[Mon(4)*U**F(2)/formula('C'),
                           Mon(1)/(Mon(4)*V**F(2)*formula('C'))]
            for k,value in enumerate(reconstructed):
                assert value==formula(common_rows[family][k],branch,j),(family,k,branch,j)
                comparisons.append(dict(family=family,endpoint=('lower','upper')[k],branch=branch,j=j))

# Substitution in one variable, including the sign required for a bound.
def sub(value,letter,replacement,direction='lower'):
    k=letters.index(letter)
    exponent=value.e[k]
    assert (exponent>=0 if direction=='lower' else exponent<=0)
    exponents=list(value.e);exponents[k]=0
    return Mon(value.c,exponents)*replacement**F(exponent)

def endpoint(family,which,branch=0,j=4):
    return formula(common_rows[family][which],branch,j)

def raw(left,right,branch=0,jleft=4,jright=4):
    return endpoint(left,1,j=jleft)/endpoint(right,0,branch=branch,j=jright)/Mon(uniform['h'])

coarse=formula('GA^7B^5S^3')
fine=formula('A^2T^5/S^2')
tgap=formula('2BS')
start=formula('(1-eps)A^2C^7/(S^2T^4)')/endpoint('Q',0)
start=sub(sub(start,'C',fine),'T',tgap)
q_n2=sub(sub(raw('Q','N_2'),'C',fine),'T',tgap)
construction=[
    start,
    q_n2,
    sub(raw('N_2','P',branch=0),'C',coarse),
    sub(raw('N_2','P',branch=1),'C',coarse),
    sub(sub(raw('P','N_3'),'C',coarse),'S',formula('B^2')),
    sub(sub(raw('P','A_3'),'C',coarse),'S',formula('B')),
    sub(sub(raw('A_3','M'),'C',coarse),'S',formula('B^2'),'upper'),
    sub(raw('M','N_3'),'C',coarse),
    sub(sub(raw('N_3','A_4'),'C',coarse),'S',formula('B^3')),
    sub(sub(raw('N_3','N_{3,bs}'),'C',coarse),'S',formula('B^3'),'upper'),
    sub(raw('N_{3,bs}','A_4'),'C',formula('rho R_{bs}A^7B^5S^3')),
    sub(sub(raw('A_4','O_j',jright=4),'C',coarse),'S',formula('B')),
    sub(sub(raw('O_j','A_5',jleft=4),'C',coarse),'S',formula('B')),
    sub(sub(raw('A_5','O_j',jright=5),'C',coarse),'S',formula('B')),
    sub(sub(raw('O_j','O_j',jleft=5,jright=6),'C',coarse),'S',formula('B')),
]
quotients=[]
for line in new_overlaps.splitlines():
    if '&' not in line or not line.lstrip().startswith((r'\(', 'START')):continue
    cells=line.rstrip().removesuffix(r'\\').split('&')
    monomial=formula(cells[1])
    idx=len(quotients)
    assert construction[idx]==monomial,(idx,construction[idx].e,monomial.e)
    assert all(e>=0 for e in monomial.e)
    q=monomial.value(dict(A=3,B=3,S=3,T=3,C=3))
    assert q>1
    quotients.append(dict(junction=cells[0].strip(),formula=cells[1].strip(),exact=str(q),decimal=float(q)))
assert len(quotients)==15
minimum=min(quotients,key=lambda z:F(z['exact']))
assert F(minimum['exact'])==F(7834457305497600,5386911336360281)>F(7,5)
assert all(F(z['exact'])>F(13,10) for z in quotients)
tail=F(400,105)*3**12
assert tail>1
assert '7834457305497600}{5386911336360281}' in candidate
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED
result=dict(status='PASS',article_sha256=EXPECTED,
            endpoint_reconstructions=comparisons,
            overlap_reconstructions=15,overlap_quotients=quotients,
            minimum=minimum,ordinary_tail_lower_bound=str(tail))
out=RESULTS/'interval_overlaps.json'
out.write_text(json.dumps(result,indent=2)+'\n')
print('PASS:',len(comparisons),'endpoint reconstructions and 15 overlap reconstructions')
print('Minimum quotient:',minimum['exact'],'>',F(7,5))
print('Infinite-tail lower bound:',tail,'> 1')
print('Output:',out.name)
