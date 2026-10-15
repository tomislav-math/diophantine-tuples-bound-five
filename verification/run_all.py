#!/usr/bin/env python3
"""Run the numerical checks accompanying the article and save their actual output."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
JOBS = [
    ('01_structural', 'Structural identities and the finite tail', 'structural_checks.py', []),
    ('02_gap_bounds', 'Explicit growth and interpolation constants', 'gap_bounds.py', []),
    ('03_approximation_polynomials', 'Approximation polynomials and coefficient tables', 'approximation_polynomials.py', []),
    ('04_approximation_bounds', 'Root-transfer and normalization estimates', 'approximation_bounds.py', []),
    ('05_small_element_bounds', 'Constants for the small-element cases', 'small_element_bounds.py', []),
    ('06_small_triples', 'The ten small-triple classes', 'extension_graphs.py', ['triples']),
    ('07_pair_graphs', 'The 48 small-pair extension graphs', 'extension_graphs.py', ['pairs']),
    ('08_pell_graphs', 'The exceptional-pair graphs in 1016 fields', 'extension_graphs.py', ['pell']),
    ('09_one_small_element', 'First pairs, third elements and cutoffs', 'one_small_element.py', []),
    ('10_interval_overlaps', 'The interval chain with no small elements', 'interval_overlaps.py', []),
    ('11_independent_intervals', 'Independent check of the remaining interval chains', 'independent_intervals.py', []),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, help='A new output directory; existing directories are refused')
    parser.add_argument('--only', choices=[x[0] for x in JOBS], help='Run one named check instead of all checks')
    parser.add_argument('--list', action='store_true', help='List the checks without running them')
    args = parser.parse_args()
    if args.list:
        for key, title, script, extra in JOBS:
            print(key + ': ' + title + ' [' + ' '.join([script, *extra]) + ']')
        return 0
    if sys.version_info < (3, 10):
        parser.error('Python 3.10 or later is required.')
    if sys.flags.optimize:
        parser.error('Do not use -O or -OO, because the checks use assertions.')
    selected = [row for row in JOBS if args.only is None or row[0] == args.only]
    source = ROOT / 'data/article.tex'
    expected = '9cb6ed2154056e158d2c420012e94ef7c397c69d66c0fd9301a1bb52197dc432'
    if digest(source) != expected:
        raise SystemExit('The article input does not match the dated paper.')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    destination = (args.results or ROOT / 'reruns' / stamp).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    environment = os.environ.copy()
    environment['D1_RESULTS_DIR'] = str(destination)
    started = datetime.now(timezone.utc).isoformat()
    records = []
    for key, title, script, extra in selected:
        program = ROOT / 'programs' / script
        before = digest(program)
        command = [sys.executable, '-u', str(program), *extra]
        print('RUN  ' + key + ' — ' + title, flush=True)
        begin = time.monotonic()
        with (destination / (key + '.txt')).open('w', encoding='utf-8') as stream:
            completed = subprocess.run(command, cwd=ROOT, env=environment,
                                       stdout=stream, stderr=subprocess.STDOUT)
        after = digest(program)
        good = completed.returncode == 0 and before == after
        records.append({'id': key, 'description': title,
                        'program': 'programs/' + script, 'arguments': extra,
                        'program_sha256': before, 'program_unchanged_during_run': before == after,
                        'output': key + '.txt', 'exit_code': completed.returncode,
                        'seconds': round(time.monotonic() - begin, 3),
                        'status': 'PASS' if good else 'FAIL'})
        print(records[-1]['status'] + ' ' + key, flush=True)
        if not good:
            break
    success = len(records) == len(selected) and all(x['status'] == 'PASS' for x in records)
    summary = {'status': 'PASS' if success else 'FAIL',
               'all_checks_run': args.only is None and len(records) == len(JOBS),
               'started_utc': started, 'finished_utc': datetime.now(timezone.utc).isoformat(),
               'python': platform.python_version(), 'platform': platform.platform(),
               'article_sha256': digest(source), 'runner_sha256': digest(Path(__file__)),
               'checks': records}
    (destination / 'run_details.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    lines = ['# Verification results', '',
             'Run started: ' + started + ' (UTC).',
             'Python: ' + platform.python_version() + '.', '',
             '| Check | Result | Program output | Time |',
             '|---|---|---|---|']
    for r in records:
        lines.append('| ' + r['description'] + ' | ' + r['status'] + ' | [Output](' + r['output'] + ') | ' + str(r['seconds']) + ' s |')
    lines += ['', ('All eleven checks passed.' if success and args.only is None else 'Overall result: ' + summary['status'] + '.'), '',
              '[Execution details and source hashes](run_details.json)', '']
    (destination / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
    print('Results: ' + str(destination), flush=True)
    return 0 if success else 1


if __name__ == '__main__':
    raise SystemExit(main())
