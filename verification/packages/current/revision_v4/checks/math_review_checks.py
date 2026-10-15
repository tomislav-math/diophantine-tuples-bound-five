#!/usr/bin/env python3
"""Read-only, hash-bound controls for the scoped v4 mathematical review.

Run from the workspace root. The article is only read. Every numerical
comparison uses integers or fractions. Human proof reconstruction is in
revision_v4/reviews/math_review.txt; these controls supplement that review.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import difflib
import json
import re

ROOT = Path(__file__).resolve().parents[2]
PACKET = ROOT / "revision_v4/checkpoints/math_packet.txt"
BASE = ROOT / "output/D1_bound5-v3.txt"
EXPECTED = "8675f3f4b0ed54cb387408844fe757fcd803f9cc01ac3ba294796b923c55881b"
source = PACKET.read_text()
base = BASE.read_text()
checks = []


def check(name, condition, **details):
    assert condition, name
    checks.append({"name": name, "pass": True,
                   **{k: str(v) for k, v in details.items()}})


def block(label):
    start = source.index("\\begin{lemma}\\label{" + label + "}")
    end = source.index("\\end{proof}", start) + len("\\end{proof}")
    return source[start:end]


def poly_mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def binomial_half(n):
    out = F(1)
    for j in range(n):
        out *= (F(1, 2) - j) / (j + 1)
    return out


check("immutable packet SHA256", sha256(PACKET.read_bytes()).hexdigest() == EXPECTED)
labels = ["lem:product-four", "lem:interpolation-parameters", "lem:direct-four",
          "lem:retained-gap", "sm:bs7-cutoff"]
blocks = {label: block(label) for label in labels}

# Bind each checked argument and its consumers to the actual source text.
needles = {
    "lem:product-four": [r"D>\frac32ABC", r"\|E\|_2<11/96",
                         r"s=59/96", r"h=155/96", r"D>6C", r"D>18",
                         r"\frac{18212679025}{27518828544}D^2<\frac23D^2"],
    "lem:interpolation-parameters": [r"(1-w)^2(1+w)=(1-w)(1-w^2)",
                                     r"(1-\rho)(1-\rho^2)\le16P^3|f_3|",
                                     r"\le(1+\rho)(1+\rho^2)",
                                     r"45/1024\le P^3|f_3|\le85/1024"],
    "lem:direct-four": [r"\frac{13BS^2T^3}{Z}+\frac{49ST^2Z}{50AP^4}",
                        r"|\Delta_4|\le24BS^2T^3", r"\frac{24\cdot101}{199}<13",
                        r"\frac{2400u}{99}<\frac{49}{50}"],
    "lem:retained-gap": [r"\frac{18}{17}A^4B^2>1", r"C>A^2T^5/S^2>16BST",
                         r"\frac{169}{16^6kA^2B^4}<\frac1{10000}",
                         r"|\sqrt{cd+1}|>(99/98)AC^4/(ST^2)"],
    "sm:bs7-cutoff": [r"\frac{39M^5}{Z}+\frac{49M^3Z}{50C^4}",
                      r"D>\frac{C^7}{100M^6}=10^{-50}C^7"],
}
for label, items in needles.items():
    for needle in items:
        check(f"source guard: {label}: {needle}", needle in blocks[label])

consumer_lines = [i for i, line in enumerate(source.splitlines(), 1)
                  if r"\ref{lem:direct-four}" in line]
check("exactly two direct-four theorem consumers", len(consumer_lines) == 2,
      lines=consumer_lines)
check("both consumers belong to reviewed lemmas",
      all(r"\ref{lem:direct-four}" in blocks[label]
          for label in ("lem:retained-gap", "sm:bs7-cutoff")))
for old in (r"\frac{9BS^2T^3}{Z}", r"\frac{4ST^2Z}{5AP^4}",
            r"\frac{27M^5}{Z}", r"\frac{4M^3Z}{5C^4}",
            r"|\Delta_4|\le16BS^2T^3", r"\frac{18}{17}A^3B^2"):
    check(f"stale altered formula absent: {old}", old not in source)

# Every textual alteration against v3 is inside one of the reviewed blocks.
source_lines = source.splitlines()
base_lines = base.splitlines()
intervals = []
for label, value in blocks.items():
    pos = source.index(value)
    begin = source[:pos].count("\n")
    intervals.append((begin, begin + value.count("\n") + 1, label))
changes = []
for tag, i, j, k, l in difflib.SequenceMatcher(None, base_lines, source_lines,
                                             autojunk=False).get_opcodes():
    if tag == "equal":
        continue
    owners = [label for begin, end, label in intervals if begin <= k and l <= end]
    check(f"diff contained in scoped block at packet line {k+1}", len(owners) == 1)
    changes.append({"kind": tag, "v3_lines": [i+1, j],
                    "packet_lines": [k+1, l], "block": owners[0]})

# Product-four: branch bounds, symmetric majorant, singular values, determinant.
check("Re square-root bounds sufficient for F entries",
      F(2, 3) > F(4, 5)**2 and F(8, 9) > F(15, 16)**2)
check("first-row entry error", F(25, 486) < F(1, 18))
check("remaining entry error", F(128, 8649) < F(1, 64))
row_a = 3 * F(1, 18) / F(3, 2)
row_other = F(3, 2) * F(1, 18) + 2 * F(1, 64)
check("majorant weighted row sums", row_a == F(1, 9) and row_other == F(11, 96))
s = F(1, 2) + F(11, 96)
h = F(3, 2) + F(11, 96)
det_coefficient = h*s*s*(1 + F(3, 36) + s/324)
check("product-four exact determinant coefficient",
      det_coefficient == F(18212679025, 27518828544), value=det_coefficient)
check("product-four implies D > (3/2)ABC", det_coefficient < F(2, 3),
      margin=F(2, 3)-det_coefficient)
check("product-five still may use older 1/8 majorant", F(11, 96) < F(1, 8))
check("cutoff's unchanged 3/8 consequence remains implied", F(3, 2) > F(3, 8))

# Three anchors: factorization and the explicit rational endpoints.
check("f3 factorization", poly_mul(poly_mul([1, -1], [1, -1]), [1, 1])
      == poly_mul([1, -1], [1, 0, -1]))
rho = F(1, 4)
check("f3 endpoint bounds", (1-rho)*(1-rho*rho)/16 == F(45, 1024)
      and (1+rho)*(1+rho*rho)/16 == F(85, 1024))
check("f3 tail prefactor", 4*abs(binomial_half(4)) == F(5, 32))
for H, lam, U in ((32, F(1, 10), F(179, 10)), (48, F(19, 200), F(15))):
    r = F(5, 32)*(F(H, H-1)**3 - 1)
    lo, hi = F(45, 1024)-r, F(85, 1024)+r
    check(f"interpolation row H={H}", lo > 0 and F(101, 199) < F(51, 100)
          and 100*hi/99 < lam and 100/(201*hi) > 5 and 100/(199*lo) < U,
          lower=100/(201*hi), upper=100/(199*lo))

# Vandermonde exponents can be sorted to (0,1,2,3) by favorable swaps.
# For ordered positive moduli, each inverted adjacent swap multiplies the
# monomial by (larger_modulus/smaller_modulus)^(positive exponent).
monomials = list(permutations(range(4)))
for exponents in monomials:
    working = list(exponents)
    for sweep in range(4):
        for j in range(3):
            if working[j] > working[j+1]:
                working[j], working[j+1] = working[j+1], working[j]
    assert working == [0, 1, 2, 3]
check("Vandermonde has 24 bounded monomials", len(monomials) == 24)
check("f4 tail prefactor", 4*abs(binomial_half(5)) == F(7, 64))
k = F(1, 1000)
r = F(7, 64)*(F(1000, 999)**4-1)
lo = (1-k)**2*(5-6*k-5*k*k)/128-r
hi = (1+k)**2*(5+6*k+5*k*k)/128+r
check("direct-four positive F4 lower bound", lo > 0, lower=lo)
check("direct-four first coefficient", F(24*101, 199) < 13,
      value=F(24*101, 199))
check("direct-four second coefficient", 2400*hi/99 < F(49, 50),
      value=2400*hi/99)
check("direct-four zero-case lower bound unchanged", 100/(201*hi) > 10,
      value=100/(201*hi))
check("direct-four zero-case upper bound unchanged", 100/(199*lo) < 17,
      value=100/(199*lo))

# Retained gap: antecedent margins, corrected A exponent and first summand.
mu = F(51051, 400000)
alpha_min = 10*(F(1, 8)-mu/9)
check("retained main coefficient exceeds one",
      F(9999, 10000)*alpha_min**2 > 1, alpha_min=alpha_min)
check("retained L3 nonvanishing margin", 18 > F(179, 10)*F(10001, 10000))
check("retained corrected A^4 nonvanishing minimum", F(18, 17)*9 > 1)
check("retained C > 16BST from T > 2ABS",
      1**6*3**3*3 >= 1, minimum_ratio=1**6*3**3*3)
k_ret = F(3999, 4000)
first_squared_max = F(169)/(16**6*k_ret*3**4)
check("retained L4 first summand < 1/100",
      first_squared_max < F(1, 10000), squared_bound=first_squared_max)
z_factor = (1-F(1, 100))/F(49, 50)
check("retained resulting Z factor", z_factor == F(99, 98))
check("retained subtract-one loss harmless",
      (z_factor**2-1)*1000**2 > 1,
      minimum_excess=(z_factor**2-1)*1000**2)

# The final branch of sm:bs7-cutoff, at its worst permitted thresholds.
M, C = 10**8, 10**100
check("cutoff P and nonvanishing hypotheses",
      C >= 4*M and C >= 1000 and C*C > 32*3*M*M)
check("cutoff preliminary first summand < 1/2", F(30*M*M, C) < F(1, 2))
check("cutoff preliminary D lower bound after subtracting one",
      (F(1, 36)-F(1, 100))*C**6/M**2 > 1)
check("cutoff preliminary D > 1000C", C**4 > 100000*M**2)
check("cutoff preliminary D > 17ABSTC^3", C**2 > 15300*M**4)
check("cutoff direct-four coefficient substitution", 13*3 == 39)
check("cutoff revised first summand < 1/2", F(78*M**5, C) < F(1, 2))
z_coefficient = F(1, 2)/F(49, 50)
check("cutoff revised Z coefficient", z_coefficient == F(25, 49))
check("cutoff revised D bound survives subtract-one",
      (z_coefficient**2-F(1, 100))*C**8/M**6 > 1,
      coefficient_margin=z_coefficient**2-F(1, 100))
check("cutoff exponent identity", 100*M**6 == 10**50)

result = {
    "status": "PASS",
    "packet": str(PACKET.relative_to(ROOT)),
    "packet_sha256": EXPECTED,
    "baseline": str(BASE.relative_to(ROOT)),
    "baseline_sha256": sha256(BASE.read_bytes()).hexdigest(),
    "scoped_blocks": {label: {"sha256": sha256(value.encode()).hexdigest(),
                               "start_line": source[:source.index(value)].count("\n")+1}
                      for label, value in blocks.items()},
    "text_changes": changes,
    "checks": checks,
    "limitations": "Scoped exact controls; human proof reconstruction is in math_review.txt. No broad theorem audit.",
}
out = ROOT / "revision_v4/checks/math_review_checks.json"
out.write_text(json.dumps(result, indent=2) + "\n")
print(f"PASS: {len(checks)} exact/source controls; immutable SHA256 {EXPECTED}")
print(out.relative_to(ROOT))
