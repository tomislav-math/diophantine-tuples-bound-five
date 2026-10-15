"""Exact focused checks for the v4 normalization and root-transfer edits.

No third-party dependencies. The determinant comparisons are finite samples,
not a substitute for the cited all-orders identity.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path, default=ROOT / 'revision_v4/D1_bound5-v4.txt')
parser.add_argument('--hash-file', type=Path, default=ROOT / 'revision_v4/release_source.sha256')
parser.add_argument('--baseline', type=Path, default=ROOT / 'output/D1_bound5-v3.txt')
args = parser.parse_args()
BASELINE = args.baseline.resolve()
EXPECTED = 'a1119af75cf623411c678ea0e9b6e7dc9a3f59a931ee67384fc98be8d8373bd5'
assert hashlib.sha256(BASELINE.read_bytes()).hexdigest() == EXPECTED
SOURCE = args.source.resolve()
SOURCE_HASH = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
RELEASE_HASH = args.hash_file.read_text().split()[0]
assert re.fullmatch('[0-9a-f]{64}', RELEASE_HASH)
assert SOURCE_HASH == RELEASE_HASH

def normalized(value):
    return re.sub(r'\s+', ' ', value).strip()

source_text = normalized(SOURCE.read_text())
required = {
    'Bennett exact location': r'''This is the unnumbered determinant identity following (2.4) in
\cite[\S2, p.~1721]{Bennett96}.''',
    'Bennett coefficient equality': r'''$\alpha_j$ and $\nu=1/2$, his coefficients satisfy $P_{ij}(1/N)=f_{ij}$,''',
    'Bennett term normalization': r'''$(1+\alpha_j/N)^{k-h_j}N^{-h_j}=N^{-k}(N+\alpha_j)^{k-h_j}$''',
    'Bennett extra pole': r'''His integral (2.1) already includes the extra pole factor
$(w-\alpha_i)^{-1}$ in \eqref{eq:ordinary-contour}.''',
    'Bennett zero node guard': r'''the identity extends to complex interpolation points with $\alpha_0=0$.''',
    'root nonzero guards': r'''Here $z\ne0$, since otherwise
$|d|=1/C<1$, and $\theta_j\ne0$, since $A_jC>1$.''',
    'root factorization': r'''u_j^2=\theta_j^2\left(1+\frac{c-a_j}{a_jz^2}\right).''',
    'root sign choice': r'''Write $u_j=\theta_jw_j$ and choose the sign so that
$\operatorname{Re}w_j\ge0$. Then $|w_j+1|\ge1$''',
    'root bound chain': r'''|\theta_j-u_j|
 =|\theta_j|\,|w_j-1|
 \le\frac{|\theta_j|\,|c-a_j|}{A_jZ^2}
 \le\frac{Ce(A_j)}{A_jZ^2}
 \le\frac{Ce(S)}{A_jZ^2}.''',
    'root proof reuse': r'''The first inequality is the root estimate in the proof of
Lemma~\ref{lem:ordinary-residue-construction}.''',
    'root basic constant': r'''\frac{1001}{1000}\left(\frac54\right)^2<\left(\frac32\right)^2''',
    'root sharper constant': r'''\frac{1001}{1000}\left(\frac{201}{200}\right)^2
 <\left(\frac{101}{100}\right)^2''',
    'root zero-root allowance': r'''Neither $r_x$ nor $v_x$ is divided out.''',
    'R2 simplified expression': r'''R_2=\frac{27B(AC-1)^2}{16A\delta^2}''',
}
for name, snippet in required.items():
    assert normalized(snippet) in source_text, name
assert normalized(r'''R_2=\frac{27AB(AC-1)^2}{16A^2\delta^2}''') not in source_text
assert normalized(r'''Indeed, the square of the quotient in
\eqref{eq:ff-error} is''') not in source_text

def binom(x, k):
    value = F(1)
    for r in range(k):
        value *= F(x-r, r+1)
    return value

def compositions(total, n):
    if n == 1:
        yield (total,)
    else:
        for h in range(total+1):
            for tail in compositions(total-h, n-1):
                yield (h,) + tail

def determinant(matrix):
    a = [row[:] for row in matrix]
    result = F(1)
    for j in range(len(a)):
        p = next((i for i in range(j, len(a)) if a[i][j]), None)
        if p is None:
            return F(0)
        if p != j:
            a[j], a[p] = a[p], a[j]
            result = -result
        pivot = a[j][j]
        result *= pivot
        for i in range(j+1, len(a)):
            ratio = a[i][j]/pivot
            for k in range(j+1, len(a)):
                a[i][k] -= ratio*a[j][k]
    return result

def residue_matrix(nodes, k, N, original):
    size = len(nodes)
    matrix = []
    for i in range(size):
        row = []
        for j in range(size):
            value = F(0)
            for hs in compositions(k+int(i == j)-1, size):
                term = binom(F(2*k+1, 2), hs[j])
                if original:
                    term *= (1+F(nodes[j], N))**(k-hs[j])*F(1, N)**hs[j]
                else:
                    term *= F(N)**(-k)*F(N+nodes[j])**(k-hs[j])
                for ell in range(size):
                    if ell != j:
                        exponent = -k-int(i == ell)-hs[ell]
                        term *= binom(-k-int(i == ell), hs[ell])*F(nodes[j]-nodes[ell])**exponent
                value += term
            row.append(value)
        matrix.append(row)
    return matrix

samples = []
for nodes, orders, N in [((0, 2, 7), range(1, 6), 101),
                         ((-2, 0, 3, 7), range(1, 4), 103),
                         ((0, -3, 5), range(1, 4), -101)]:
    m = len(nodes)-1
    for k in orders:
        original = residue_matrix(nodes, k, N, True)
        normalized = residue_matrix(nodes, k, N, False)
        assert original == normalized
        actual = determinant(normalized)
        expected = F(N)**(-(m+1)*k)
        for v in range(-1, m):
            expected *= binom(F(1, 2)-v*k, k)
        for ell in range(m+1):
            for r in range(m+1):
                if r != ell:
                    expected *= F(nodes[r]-nodes[ell])**(-k)
        assert actual == expected != 0
        samples.append({'nodes': nodes, 'k': k, 'N': N,
                        'determinant': str(actual), 'normalization_matches': True})

# Clear x^2*c*(c*d+1). Polynomial coefficients use monomials (x,c,d).
def add(*args):
    out = {}
    for p in args:
        for mon, v in p.items():
            out[mon] = out.get(mon, 0)+v
    return {mon: v for mon, v in out.items() if v}

def mul(p, q):
    out = {}
    for u, a in p.items():
        for v, b in q.items():
            mon = tuple(i+j for i, j in zip(u, v))
            out[mon] = out.get(mon, 0)+a*b
    return out

x, c, d, one = ({(1,0,0):1}, {(0,1,0):1}, {(0,0,1):1}, {(0,0,0):1})
minus_x = {(1,0,0):-1}
left = mul(c, mul(add(mul(x,c), one), add(mul(x,d), one)))
right = mul(add(mul(x,c), one), add(mul(x, add(mul(c,d), one)), c, minus_x))
assert left == right
assert F(1001, 1000)*F(5, 4)**2 < F(3, 2)**2
assert F(1001, 1000)*F(201, 200)**2 < F(101, 100)**2

# This tests the cancellation independently at exact rational values.
for A, B, C, delta in [(F(1), F(3), F(1001), F(2)),
                        (F(7, 2), F(19, 4), F(900), F(11, 3))]:
    assert F(27)*A*B*(A*C-1)**2/(16*A*A*delta*delta) == F(27)*B*(A*C-1)**2/(16*A*delta*delta)

result = {
    'status': 'PASS', 'baseline_sha256': EXPECTED,
    'source': str(SOURCE), 'source_sha256': SOURCE_HASH,
    'release_hash_matches': True,
    'source_controls': {name: True for name in required},
    'old_R2_expression_removed': True,
    'duplicated_root_proof_removed': True,
    'Bennett_source': 'https://personal.math.ubc.ca/~bennett/B-Binom.pdf',
    'Bennett_location': 'Section 2, p. 1721, unnumbered determinant identity following (2.4)',
    'Bennett_samples': samples,
    'root_transfer_cleared_polynomial_identity': True,
    'root_transfer_constants_squared': {
        'basic': str(F(1001,1000)*F(5,4)**2),
        'sharper': str(F(1001,1000)*F(201,200)**2)},
    'R2_cancellation': True,
    'scope': 'Exact finite normalization/determinant checks and polynomial identity. The cited theorem supplies the all-orders determinant. Analytic sign and nonzero arguments were reviewed in the report.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
print('PASS: release source/hash and 16 text controls; 11 Bennett normalization/determinant samples, root identity, both root constants, and R2 cancellation.')
