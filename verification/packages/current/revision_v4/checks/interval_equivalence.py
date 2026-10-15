#!/usr/bin/env python3
"""Exact, source-based checks of the scoped interval/M revision.

No manuscript is modified. The default candidate is the noncanonical
snapshot. Pass a candidate filename to check root's assembled revision.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'output/D1_bound5-v3.txt'
CANDIDATE = Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'revision_v4/proposals/interval_proposed_snapshot.txt'
base, candidate = BASE.read_text(), CANDIDATE.read_text()
assert hashlib.sha256(BASE.read_bytes()).hexdigest() == 'a1119af75cf623411c678ea0e9b6e7dc9a3f59a931ee67384fc98be8d8373bd5'

def table_at(text, marker, before=False):
    pos=text.index(marker)
    begin=text.rfind(r'\begin{tabular}',0,pos) if before else text.index(r'\begin{tabular}',pos)
    end=text.index(r'\end{tabular}',begin)+len(r'\end{tabular}')
    return text[begin:end]

old_shared=table_at(base,r'\label{tab:pr-windows}')
new_shared=table_at(candidate,r'\label{tab:pr-windows}')
assert old_shared==new_shared
old_overlaps=table_at(base,r'\label{tab:pr-overlaps}')
new_overlaps=table_at(candidate,r'\label{tab:pr-overlaps}')
assert old_overlaps==new_overlaps
old6=table_at(base,r'\paragraph{The approximation intervals.}')
old9=table_at(base,'For completeness, all raw interval endpoints used here are displayed')
assert 'tab:fixed56-windows' not in candidate
assert candidate.count(r'\label{tab:pr-windows}')==1
assert candidate.count(r'\label{tab:pr-overlaps}')==1
assert candidate.index(r'\label{tab:pr-windows}') < candidate.index(r'\section{Entries of modulus at least three}')
assert candidate.index(r'\label{tab:pr-overlaps}') < candidate.index(r'\section{Entries of modulus at least three}')

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

common_rows, section6_rows, section9_rows = rows(new_shared),rows(old6),rows(old9)
comparisons=[]
for key, pair in common_rows.items():
    for k,s in enumerate(pair):
        for branch in (0,1) if key=='P' else (0,):
            for j in (4,5,6) if key=='O_j' else (4,):
                six=section6_rows[key][k].replace('735','243')
                assert formula(s,branch,j)==formula(six,branch,j),(key,k,'section6')
                assert formula(s,branch,j,unit_A=(k==1))==formula(section9_rows[key][k],branch,j),(key,k,'section9')
                comparisons.append(dict(row=key,endpoint=('l','u')[k],branch=branch,j=j))

# These endpoint comparisons are formal Laurent monomial equalities, not
# evaluations at sample moduli. They also compare both branches of max(P).
# The ordinary formula is identical in j; 4,5,6 additionally test its use.

quotients=[]
for line in new_overlaps.splitlines():
    if '&' not in line or not line.lstrip().startswith((r'\(', 'START')):continue
    cells=line.rstrip().removesuffix(r'\\').split('&')
    q=formula(cells[1]).value(dict(A=3,B=3,S=3,T=3,C=3))
    assert q>1
    quotients.append(dict(junction=cells[0].strip(),exact=str(q),decimal=float(q)))
assert len(quotients)==15
minimum=min(quotients,key=lambda z:F(z['exact']))
assert F(minimum['exact'])==F(7834457305497600,5386911336360281)>F(7,5)
assert all(F(z['exact'])>F(13,10) for z in quotients)
tail=F(400,105)*3**12
assert tail>1
assert '7834457305497600}{5386911336360281}' in candidate

# The specialized 13-cutoff table and its N0 prescription remain verbatim.
cut_marker='The next table gives a strict lower bound $c_ex^{p_e}$'
old_cut=base[base.index(cut_marker):base.index('If $N(t)\\ge N_0$')]
new_cut=candidate[candidate.index(cut_marker):candidate.index('If $N(t)\\ge N_0$')]
assert old_cut==new_cut
graph='7 & $(-5,-1)$ & $(0,2)$ & 737 & 735'
assert graph in base and graph in candidate
for lab in ('tab:pr-forty-coefficients','tab:pr-five-coefficients','tab:fixed56'):
    # These tables have captions/labels before OR after their tabular block.
    pos=base.index('\\label{'+lab+'}')
    env=max(('table','longtable'),key=lambda e:base.rfind('\\begin{'+e+'}',0,pos))
    prev=base.rfind('\\begin{'+env+'}',0,pos)
    ending='\\end{'+env+'}'
    end=base.index(ending,pos)+len(ending)
    table=base[prev:end]
    assert table in candidate,lab

result=dict(base_sha256=hashlib.sha256(BASE.read_bytes()).hexdigest(),
            candidate=str(CANDIDATE),candidate_sha256=hashlib.sha256(CANDIDATE.read_bytes()).hexdigest(),
            common_endpoint_table_unchanged=True,common_overlap_table_unchanged=True,
            formal_endpoint_comparisons=comparisons,section6_source_parsed_quotients=quotients,
            minimum=minimum,ordinary_tail_lower_bound=str(tail),
            specialized_cutoff_formula_and_tables_unchanged=True,unrelated_735_count_preserved=True,
            all_checks_passed=True)
out=Path(__file__).with_suffix('.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_checks_passed':True,'formal_endpoint_comparisons':len(comparisons),'quotients':len(quotients),'minimum':minimum,'output':str(out)},indent=2))
