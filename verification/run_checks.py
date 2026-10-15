#!/usr/bin/env python3
"""Reproduce the article's numerical checks in a temporary working copy."""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--quick", action="store_true", help="Run seven scoped exact checks; omit large graphs and full polynomial expansion")
    mode.add_argument("--full", action="store_true", help="Also regenerate the polynomial families and all graph/closure computations")
    parser.add_argument("--results", type=Path, help="New directory for logs and regenerated results")
    args = parser.parse_args()
    if sys.version_info < (3, 10):
        parser.error("Python 3.10 or later is required; no third-party packages are needed.")
    if sys.flags.optimize:
        parser.error("Run without -O or -OO: the checks use assertions.")

    manifest = json.loads((ROOT / "MANIFEST.json").read_text())
    for name, expected in manifest["sha256"].items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit("Integrity check failed: " + name)
    print(f"PASS: {len(manifest['sha256'])} file hashes", flush=True)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    destination = (args.results or ROOT / "runs" / stamp).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    current = "current"
    enumeration = "enumeration"
    jobs = [
        ("structural", enumeration, ["revision_v2/checks/structural_v2.py", "--expected-sha", "5edcf8777371c4ad8979f122022d3119a8ae2de75d9674b336df994bba860b70"]),
        ("small_element_constants", enumeration, ["revision_v2/checks/small_v2_constants.py"]),
        ("current_structural_constants", current, ["revision_v4/checks/math_review_checks.py"]),
        ("root_transfer_and_normalization", current, ["revision_v4/checks/approximation_local.py"]),
        ("one_small_element", current, ["revision_v4/checks/completion_v4.py"]),
        ("interval_equivalence", current, ["revision_v4/checks/interval_equivalence.py", "revision_v4/D1_bound5-v4.txt"]),
        ("independent_intervals", current, ["revision_v4/reviews/interval_referee_independent.py"]),
    ]
    if args.full:
        jobs += [("polynomial_families", enumeration, ["revision_v2/checks/approximation_exact.py"])]
        jobs += [("small_" + name, enumeration, ["revision_v2/checks/small_v2_exact.py", name])
                 for name in ("triples", "pairs", "pell")]
    results = []
    with tempfile.TemporaryDirectory(prefix="diophantine_tuples_checks_") as temporary:
        work = Path(temporary)
        shutil.copytree(ROOT / "packages", work / "packages")
        for name, package, command in jobs:
            print("RUN " + name, flush=True)
            start = time.monotonic()
            with (destination / (name + ".log")).open("w") as output:
                completed = subprocess.run([sys.executable, "-u", *command], cwd=work / "packages" / package,
                                           stdout=output, stderr=subprocess.STDOUT)
            record = {"check": name, "package": package, "command": command,
                      "exit_code": completed.returncode, "seconds": round(time.monotonic() - start, 3)}
            results.append(record)
            print(("PASS " if completed.returncode == 0 else "FAIL ") + name, flush=True)
            if completed.returncode:
                print("See " + str(destination / (name + ".log")), file=sys.stderr)
                break
        for item in (work / "packages").rglob("*.json"):
            target = destination / "generated" / item.relative_to(work / "packages")
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, target)
    success = len(results) == len(jobs) and all(item["exit_code"] == 0 for item in results)
    summary = {"status": "PASS" if success else "FAIL", "mode": "full" if args.full else "quick",
               "python": sys.version, "checks": results,
               "large_graphs_recomputed": bool(args.full and success),
               "scope": "Numerical and symbolic checks listed in README.md; a successful run does not replace the analytic proofs in the article."}
    (destination / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("Results: " + str(destination), flush=True)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
