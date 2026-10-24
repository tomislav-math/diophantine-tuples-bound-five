# Numerical verification of the article

This directory accompanies [the paper dated 15 October 2026](../Diophantine_tuples_bound_five_2026-10-15.pdf).
The table below identifies each calculation, the method used, its program,
and the actual output of the program. The mathematical checks use exact integer,
rational or polynomial arithmetic; decimal values in the output are for display.

**All eleven checks passed in a complete run on 15 October 2026.**
The full run takes about three minutes on the recorded machine; timings
for each calculation are given on the results page.

## What is checked

| Part of the paper | What the program does | Program | Output |
|---|---|---|---|
| Lemmas 2.2 and 2.5 | Enumerates the finite tail candidates and checks the structural identities and determinant constants | [structural_checks.py](programs/structural_checks.py) | [01_structural.txt](results/01_structural.txt) |
| Explicit bounds in Sections 2–3 and Lemma 7.4 | Checks the rational inequalities used in the growth, interpolation and fourth-element estimates, including those used in the two-small-entry cutoff $C<10^{550}$ | [gap_bounds.py](programs/gap_bounds.py) | [02_gap_bounds.txt](results/02_gap_bounds.txt) |
| Approximation constructions in Sections 4–5 | Constructs the polynomial matrices and verifies their coefficients, divisibility, determinants and residual cancellations | [approximation_polynomials.py](programs/approximation_polynomials.py) | [03_approximation_polynomials.txt](results/03_approximation_polynomials.txt) |
| Normalization and root transfer in Section 4 | Checks the root-transfer algebra and bounds, together with explicitly listed finite determinant examples | [approximation_bounds.py](programs/approximation_bounds.py) | [04_approximation_bounds.txt](results/04_approximation_bounds.txt) |
| Constants in Section 7 | Verifies the numerical inequalities used to obtain the finite search bounds | [small_element_bounds.py](programs/small_element_bounds.py) | [05_small_element_bounds.txt](results/05_small_element_bounds.txt) |
| Lemma 7.12 and Table 5 | Enumerates the small triples and computes their extension closures through modulus 10^40 | [extension_graphs.py triples](programs/extension_graphs.py) | [06_small_triples.txt](results/06_small_triples.txt) |
| Lemma 7.6 and Table 4 | Classifies the small pairs, constructs all 48 extension graphs through modulus 10^550, and tests every graph for triangles | [extension_graphs.py pairs](programs/extension_graphs.py) | [07_pair_graphs.txt](results/07_pair_graphs.txt) |
| Lemma 7.9 | Uses continued fractions to enumerate the extensions of the pair {-1,1} in 1,016 fields and checks the resulting graphs | [extension_graphs.py pell](programs/extension_graphs.py) | [08_pell_graphs.txt](results/08_pell_graphs.txt) |
| Sections 8–9, Tables 6–10 | Enumerates the first pairs and third elements, regenerates the specialized coefficients, and verifies the interval overlaps and final cutoffs | [one_small_element.py](programs/one_small_element.py) | [09_one_small_element.txt](results/09_one_small_element.txt) |
| The interval chain in Sections 5.5–6 | Checks the endpoint formulas and the overlaps used when no entry is small | [interval_overlaps.py](programs/interval_overlaps.py) | [10_interval_overlaps.txt](results/10_interval_overlaps.txt) |
| Interval calculations in Sections 8–9 | Reconstructs the coefficients, overlaps and fixed-triple cutoffs with a separate implementation | [independent_intervals.py](programs/independent_intervals.py) | [11_independent_intervals.txt](results/11_independent_intervals.txt) |

The [results page](results/README.md) gives the status and running time of
every check. Detailed numerical results are saved as JSON files beside the
text outputs. The complete execution record, including the hashes of the
programs actually run, is in [run_details.json](results/run_details.json).

In the graph outputs, a coordinate pair `(u, v)` represents
$(u+v\sqrt{-d})/2$. In the one-small-entry calculations, it represents
$u+v\omega$, with $\omega$ as defined in the article. Junction numbers
start at 1.

## How to reproduce the results

Install **Python 3.10 or later**. No additional Python packages are required.
Download the repository with **Code → Download ZIP**, then extract it.
Open a terminal in the extracted repository folder and run:

```sh
python3 verification/run_all.py
```

This runs all eleven checks, including all large graph computations.
It creates a new directory under `verification/reruns/` and prints its
location. The supplied results are not overwritten. Each program's complete
output and exit status are saved. The command succeeds only if every check
succeeds.

To run just one calculation, for example the 48 pair graphs:

```sh
python3 verification/run_all.py --only 07_pair_graphs
```

Use `--list` to list the check names, or `--results PATH` to select a new
output directory. Run Python without `-O` or `-OO`, since the programs use
assertions to test the identities and inequalities.

## Input and interpretation

The only manuscript input is [data/article.tex](data/article.tex), the complete
LaTeX source of the PDF linked above. It supplies the printed formulas and
tables against which the calculations are compared. The programs perform
the calculations themselves; saved result files are not used as answers.
The source is read as text. No LaTeX installation or typesetting software is required.

The programs verify the finite computations and the algebraic identities
listed above. The proofs that the bounded searches cover all relevant
possibilities are given in the article. The ordinary-family determinant examples
are finite tests; the general identity is the cited theorem used in the paper.
