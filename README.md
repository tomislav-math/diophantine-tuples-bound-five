# A uniform bound of five for Diophantine tuples in imaginary quadratic integer rings

**Zrinka Franušić and Tomislav Pejković**

15 October 2026

**[Read the paper (PDF)](Diophantine_tuples_bound_five_2026-10-15.pdf)**

A $D(1)$-tuple is a finite set of distinct nonzero elements such that the
product of any two distinct elements, increased by one, is a square.
We prove that every $D(1)$-tuple in the full ring of integers of an
imaginary quadratic field has at most five elements. Zero is allowed
as a square in this definition.

The proof combines determinant estimates, interpolation of square roots
and explicit simultaneous approximation. The remaining cases are reduced
to bounded computations in exact integer and rational arithmetic. The
article gives the search bounds and the arguments establishing that the
enumerations are complete.

The bound leaves open whether a $D(1)$-quintuple exists in the ring of
integers of an imaginary quadratic field.

## Files

| File or directory | Contents |
|---|---|
| [Diophantine_tuples_bound_five_2026-10-15.pdf](Diophantine_tuples_bound_five_2026-10-15.pdf) | The article, 58 pages |
| [verification/](verification/) | Numerical verification programs, input data and recorded results |
| [verification/README.md](verification/README.md) | What each program checks, requirements and reproduction commands |

The verification guide lists each calculation, explains its method, and
links directly to the program and its recorded output. All eleven groups
of checks passed in a complete run on 15 October 2026. The full verification
can be repeated with one command and requires only Python.

The code in verification/ is released under the MIT License; the article (PDF and LaTeX source) is not covered by this license and remains copyright the authors.

## Lean formalization

The main theorem has been fully formalized in Lean 4 with mathlib. The
[Lean source code and verification instructions](https://github.com/tomislav-math/d1-bound5-lean)
are available in a separate repository. The
[Prove2Me project](https://prove2.me/missions/458763a2-e999-4896-81fd-7915d24557a0)
contains the verified proof, and its
[main theorem](https://prove2.me/theorems/4bb5d407-713b-46ac-967d-97d7abf88598)
has status **Proved**. Verification on Prove2Me was completed on 24 October 2026.
