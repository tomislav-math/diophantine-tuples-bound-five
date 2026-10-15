# Numerical verification

This supplement contains the exact-arithmetic programs used to check the numerical and symbolic computations in **A uniform bound of five for Diophantine tuples in imaginary quadratic integer rings**. It includes reproducible computations, their saved outputs, and earlier numerical implementations.

## Running the checks

Use **Python 3.10 or later**. The main checks use only the Python standard library. Run these commands from this directory, without `-O` or `-OO`.

```sh
python3 run_checks.py --quick
python3 run_checks.py --full
```

`--quick` runs seven exact checks, including the small first-pair and triple enumerations, coefficients, cutoffs, and interval comparisons. It omits the large extension graphs and full symbolic expansion of the approximation families. `--full` adds that expansion and all three graph/closure computations below. The full run can take substantially longer.

Each run first checks the supplied file hashes, then executes the programs in a temporary copy. Saved evidence is never overwritten. Fresh logs and generated JSON files are written under `runs/`; choose another new directory with `--results PATH`. A failed check gives a nonzero exit status. `run_summary.json` records exactly which checks ran and their exit statuses.

## Programs and coverage

| Computation | Program | Run |
| --- | --- | --- |
| Structural identities, finite irregular-tail computation (90 pairs, 118 retained triples, 16 candidates), and determinant constants | [`structural_v2.py`](packages/enumeration/revision_v2/checks/structural_v2.py) | quick/full |
| Rational bounds and identities for cases with small elements | [`small_v2_constants.py`](packages/enumeration/revision_v2/checks/small_v2_constants.py) | quick/full |
| Revised structural and interpolation constants, product bound and fourth-element estimates | [`math_review_checks.py`](packages/current/revision_v4/checks/math_review_checks.py) | quick/full |
| Root-transfer identities and normalization; finite sample checks of the determinant identity | [`approximation_local.py`](packages/current/revision_v4/checks/approximation_local.py) | quick/full |
| Forty first-pair classes, 57 triples, specialized coefficients, 41 overlap certificates, 56 cutoffs and the triple (1,3,8) | [`completion_v4.py`](packages/current/revision_v4/checks/completion_v4.py) | quick/full |
| Equivalence of the printed interval formulas and the unified interval table | [`interval_equivalence.py`](packages/current/revision_v4/checks/interval_equivalence.py) | quick/full |
| Independent reconstruction of coefficients, interval overlaps and fixed-triple cutoffs | [`interval_referee_independent.py`](packages/current/revision_v4/reviews/interval_referee_independent.py) | quick/full |
| Five finite approximation families, determinants, coefficient sums, residuals, and the additional P and Q rows | [`approximation_exact.py`](packages/enumeration/revision_v2/checks/approximation_exact.py) | full |
| Ten small-triple classes and their extension closures through modulus 10^40 | [`small_v2_exact.py triples`](packages/enumeration/revision_v2/checks/small_v2_exact.py) | full |
| Seventy-two small-pair classes and all 48 nondegenerate extension graphs through modulus 10^550 | [`small_v2_exact.py pairs`](packages/enumeration/revision_v2/checks/small_v2_exact.py) | full |
| Pell enumeration for 1,016 fields and the graphs for the pair {-1,1}, through modulus 10^550 | [`small_v2_exact.py pell`](packages/enumeration/revision_v2/checks/small_v2_exact.py) | full |

Expected graph totals are 36,931 vertices, 36,854 edges and no triangles for the 48 pair graphs. The Pell computation gives 409 nonempty graphs, 19,612 vertices, 9,806 edges and no triangles. All acceptance decisions in the main programs use exact integers or rational numbers; decimal diagnostics are for display.

## Files and provenance

- [`packages/current/`](packages/current/) contains checks of the final mathematical version.
- [`packages/enumeration/`](packages/enumeration/) contains the original polynomial and enumeration programs, retaining their internal version names and frozen inputs. Later changes to constants and interval formulas are checked by the current programs.
- [`recorded_results/`](recorded_results/) contains outputs saved during the preparation and verification of the article. Their presence does not mean that a new run has recomputed them.
- [`archive/numerical_history.zip`](archive/numerical_history.zip) preserves earlier numerical implementations and their available outputs, with a separate guide. These are historical checks and are not part of the two commands above.
- [`PROVENANCE.json`](PROVENANCE.json) records the original and distributed hashes. [`MANIFEST.json`](MANIFEST.json) records the hashes checked by the runner.

Some programs extract printed tables and formulas from frozen text snapshots. The `.txt` files inside `packages/` preserve the original input bytes, including mathematical markup; they are data inputs and require no typesetting software. Only their filenames and the corresponding pathname strings in the programs have changed. This preserves the programs' original source-hash checks.

The seven `--quick` checks passed again when this supplement was packaged on 15 October 2026; the fresh logs are in [`recorded_results/package_check_2026-10-15/`](recorded_results/package_check_2026-10-15/). The large graphs and full polynomial expansion were not rerun during this packaging step.

The scripts check the finite computations and explicitly listed algebraic identities. The completeness of the bounded enumerations uses the descent and continued-fraction arguments proved in the article. Finite sample determinant tests do not prove an identity for arbitrary order.
