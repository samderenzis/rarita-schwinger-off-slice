#!/usr/bin/env python3
"""Owner re-run of the rarita-schwinger-off-slice scripts (phase1 snapshot, 2026-11-22 22:00).

What it does
  1. Copies the off-slice repository into a fresh work folder (your copy is never modified).
  2. Checks every file against the repository's own SHA256SUMS.txt.
  3. Re-runs every script that has a committed log, in the copy, writing NEW logs and NEW data files.
  4. Compares each new log with the committed log after removing lines that are expected to differ
     (date, Python version, the command line, timing in brackets).
  5. Compares each regenerated data file (JSON) with the committed one, byte for byte (SHA-256).
  6. Writes SUMMARY.txt and MANIFEST.json, and exits 0 only if every job is MATCH.

Statuses per job
  MATCH          new output equals the committed output (after normalisation)
  MATCH-EOL      the same, except regenerated JSON has Windows line endings (CRLF); content identical
  MATCH-KNOWN    only a known, explained difference (listed in KNOWN_DIFFS)
  DIFFER         the script ran, but its output differs: inspect the diff in the job folder
  FAIL           the script exited non-zero or timed out
  NOT RUN        a prerequisite is missing (for example --parent was not given)

A full pass shows the scripts are deterministic and portable to this machine. It does not show the
mathematics is correct, and it is not independent review.

Usage (from the package folder):
  python -I offslice_rerun.py --repo <path to rarita-schwinger-off-slice> [--parent <path to rarita-schwinger-well-posedness v1.0.11>] [--skip-long]
"""
from __future__ import annotations

import argparse
import datetime as _dt
import difflib
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time

PKG_VERSION = 'offslice-owner-rerun 1.1 (2026-11-25)'
EXPECTED_SUMS_SHA256 = None  # filled from the repo; reported, not enforced

POINTS = [('A', 'A'), ('A-', 'Am'), ('B', 'B'), ('B-', 'Bm')]

# Differences found when this package was tested (2026-11-25) and explained. A job whose ONLY differences are
# exactly these line pairs (old committed line -> new line) is reported as MATCH-KNOWN, never silently as MATCH.
KNOWN_DIFFS = {
    'M2 extended certificate': {
        'why': ('the committed log was written before m2_box_tools.py gained the ray-2 / end-box tools '
                '(CHANGELOG, endpoint-gap entry); the certificate it writes is byte-identical (sha256 7bcf43c6...)'),
        'pairs': [('b5ac2af05526dae644bf2bf3ca17c9b0d8e3e0e38cc7867d4593a2370dfe5cd9  m2_box_tools.py',
                   '4c6ed9b6568a5bcbe0ee07c3ba247e9b34c78dbb6fb6a78a688b4f99849d1210  m2_box_tools.py')],
    },
}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------
# Normalisation: lines expected to differ between machines/runs
# ------------------------------------------------------------------------------------------------
_DROP_PREFIXES = ('date (utc)', 'date:', 'python ', 'command:', 'cmd:', 'run time', 'elapsed', 'wall time',
                  'host', 'platform')
_TIMING = re.compile(r'\s*[\[(]\s*[0-9]+(?:\.[0-9]+)?\s*(?:s|ms|min)\s*[\])]')
_TIMING2 = re.compile(r'\b(?:in|took|time)\s*[0-9]+(?:\.[0-9]+)?\s*(?:s|ms|min|seconds)\b', re.I)
# a token is treated as a file path only if it names a file with a known extension, or starts like a path;
# plain fractions such as 81608/27765334551 are left untouched
_PATHTOK = re.compile(r'(?:[A-Za-z]:[\\/]|\.{1,2}[\\/]|[\\/])?(?:[^\s\\/]+[\\/])+[^\s\\/]*\.(?:py|json|txt|md|cff)\b'
                      r'|(?:[A-Za-z]:|\.{1,2})?[\\/][^\s]+[\\/][^\s]*')


def normalise(text: str) -> list:
    out = []
    for line in text.replace('\r\n', '\n').split('\n'):
        s = line.strip()
        if not s:
            continue
        low = s.lower()
        if low.startswith(_DROP_PREFIXES):
            continue
        s = _TIMING.sub('', s)
        s = _TIMING2.sub('<time>', s)
        s = _PATHTOK.sub(lambda m: re.split(r'[\\/]', m.group(0))[-1], s).strip()
        if s:
            out.append(s)
    return out


# ------------------------------------------------------------------------------------------------
# Jobs
# ------------------------------------------------------------------------------------------------

def build_jobs(parent: str | None, skip_long: bool) -> list:
    J = []

    def job(name, argv, log=None, regen=None, needs=None, long=False, alt_logs=()):
        J.append({'name': name, 'argv': argv, 'log': log, 'regen': regen or [], 'needs': needs,
                  'long': long, 'alt_logs': list(alt_logs)})

    job('points self-check', ['scripts/points.py'])
    job('witness check (rebuilt parent witness_min on rs_model)', ['scripts/witness_check.py'])
    job('step 1 slice regression vs parent v1.0.11',
        ['scripts/step1_slice_regression.py', '--parent', parent or '<missing>', '--out', '@LOG'],
        log='step1_slice_regression_2026-11-22.txt', needs='parent')
    job('step 2 reduction with q free', ['scripts/step2_reduction_q_free.py', '--out', '@LOG'],
        log='step2_reduction_q_free_2026-11-22.txt')
    job('step 4 dispersion and crossings', ['scripts/step4_dispersion_and_crossings.py', '--out', '@LOG'],
        log='step4_dispersion_and_crossings_2026-11-22.txt')
    for p, tag in POINTS:
        job(f'step 5 exact point {p}',
            ['scripts/step5_point_exact.py', '--point', p, '--out', '@LOG', '--json', f'@REGEN:point_{tag}_chain.json'],
            log=f'step5_point_{tag}_exact_2026-11-22.txt', regen=[f'point_{tag}_chain.json'])
        job(f'step 5 curved symbol point {p}', ['scripts/step5_curved_symbol.py', '--point', p, '--out', '@LOG'],
            log=f'step5_point_{tag}_curved_symbol_2026-11-22.txt')
        job(f'second engine point {p}', ['scripts/verify_point_fractions.py', f'data/point_{tag}_chain.json', '--out', '@LOG'],
            log=f'verify_point_{tag}_fractions_2026-11-22.txt')
    job('M2 rotation invariance', ['scripts/m2_rotation_invariance.py', '--out', '@LOG'],
        log='m2_rotation_invariance_2026-11-22.txt')
    job('M2 q-direction certificate', ['scripts/m2_family_q_direction.py', '--out', '@LOG', '--json', '@REGEN:m2_minor.json'],
        log='m2_family_q_direction_2026-11-22.txt', regen=['m2_minor.json'])
    job('M2 second engine', ['scripts/verify_m2_fractions.py', 'data/m2_minor.json', '--out', '@LOG'],
        log='verify_m2_fractions_2026-11-22.txt')
    job('M2 extended certificate', ['scripts/m2_family_extended.py', '--out', '@LOG', '--json', '@REGEN:m2_extended_certificate.json'],
        log='m2_family_extended_2026-11-22.txt', regen=['m2_extended_certificate.json'])
    job('M2 extended second check', ['scripts/verify_m2_extended_fractions.py', 'data/m2_extended_certificate.json', '--out', '@LOG'],
        log='verify_m2_extended_fractions_2026-11-22.txt',
        alt_logs=['regression_verify_m2_extended_new_boxtools_2026-11-22.txt'])
    job('M2 endpoint gaps', ['scripts/m2_endpoint_gaps.py', '--out', '@LOG', '--json', '@REGEN:m2_gap_certificate.json'],
        log='m2_endpoint_gaps_2026-11-22.txt', regen=['m2_gap_certificate.json'])
    job('M2 endpoint gaps second check', ['scripts/verify_m2_gaps_fractions.py', 'data/m2_gap_certificate.json', '--out', '@LOG'],
        log='verify_m2_gaps_fractions_2026-11-22.txt')
    job('M2 dispersion export', ['scripts/m2_dispersion_export.py', '--out', '@LOG', '--json', '@REGEN:m2_dispersion.json'],
        log='m2_dispersion_export_2026-11-22.txt', regen=['m2_dispersion.json'])
    job('M2 dispersion fractions proof (long)', ['scripts/verify_m2_dispersion_fractions.py', 'data/m2_dispersion.json', '--out', '@LOG'],
        log='verify_m2_dispersion_fractions_2026-11-22.txt', long=True)
    job('M2 curved symbol (general)', ['scripts/m2_curved_symbol_general.py', '--out', '@LOG'],
        log='m2_curved_symbol_general_2026-11-22.txt')
    job('audit: coercivity constants', ['audit/coercivity_constants.py', '--out', '@LOG'],
        log='audit_coercivity_constants_2026-11-22.txt')
    job('audit: endpoint multiplicity', ['audit/endpoint_multiplicity.py', '--out', '@LOG'],
        log='audit_endpoint_multiplicity_2026-11-22.txt')
    job('audit: Jordan classification', ['audit/jordan_classification.py', '--out', '@LOG'],
        log='audit_jordan_classification_2026-11-22.txt')
    job('audit: restricted symbol (fractions)', ['audit/restricted_symbol_fractions.py', '--out', '@LOG'],
        log='audit_restricted_symbol_fractions_2026-11-22.txt')
    for j in J:
        if j['long'] and skip_long:
            j['skip'] = 'skipped with --skip-long'
        if j['needs'] == 'parent' and not parent:
            j['skip'] = 'no --parent given'
    return J


def data_status(newp: str, oldp: str) -> tuple:
    """'SAME' (byte-identical), 'EOL' (identical after converting CRLF to LF; Windows text-mode line endings),
    or 'DIFF'. Returns (status, new_sha, committed_sha, new_sha_after_eol_fix)."""
    nb, ob = open(newp, 'rb').read(), open(oldp, 'rb').read()
    hn, ho = hashlib.sha256(nb).hexdigest(), hashlib.sha256(ob).hexdigest()
    if nb == ob:
        return 'SAME', hn, ho, hn
    fixed = nb.replace(b'\r\n', b'\n')
    hf = hashlib.sha256(fixed).hexdigest()
    return ('EOL' if fixed == ob else 'DIFF'), hn, ho, hf


def compare_logs(new_path: str, committed_paths: list, known=None, hash_map=None) -> tuple:
    new = normalise(open(new_path, encoding='utf-8', errors='replace').read())
    if hash_map:
        new = [ln for ln in new]
        for a, b in hash_map.items():
            new = [ln.replace(a, b) for ln in new]
    if known:
        swap = {b: a for a, b in known['pairs']}
        new_k = [swap.get(x, x) for x in new]
        for cp in committed_paths:
            if os.path.exists(cp) and new_k == normalise(open(cp, encoding='utf-8', errors='replace').read()) and new_k != new:
                return 'KNOWN', os.path.basename(cp), []
    best = None
    for cp in committed_paths:
        if not os.path.exists(cp):
            continue
        old = normalise(open(cp, encoding='utf-8', errors='replace').read())
        if new == old:
            return True, os.path.basename(cp), []
        diff = list(difflib.unified_diff(old, new, 'committed/' + os.path.basename(cp), 'new', lineterm='', n=1))
        if best is None or len(diff) < len(best[1]):
            best = (os.path.basename(cp), diff)
    if best is None:
        return False, None, ['committed log not found']
    return False, best[0], best[1]


def recompare(repo: str, outdir: str) -> int:
    """Re-compare an earlier run's new_logs/ and regenerated_data/ with the repository, without re-running."""
    man_path = os.path.join(outdir, 'MANIFEST.json')
    if not os.path.isfile(man_path):
        print(f'NOT VERIFIED: no MANIFEST.json in {outdir}')
        return 2
    man = json.load(open(man_path, encoding='utf-8'))
    jobs = {j['name']: j for j in build_jobs('x', False)}
    lines, all_ok = [], man.get('sha256sums_ok', False)
    for rec in man['jobs']:
        j = jobs.get(rec['job'])
        st = rec['status']
        if st in ('MATCH', 'MATCH-KNOWN', 'MATCH-EOL', 'DIFFER') and j is not None:
            st, hash_map, eol = 'MATCH', {}, False
            for r in j['regen']:
                newp, oldp = os.path.join(outdir, 'regenerated_data', r), os.path.join(repo, 'data', r)
                if not os.path.exists(newp):
                    st = 'DIFFER'
                    continue
                ds, hn, ho, hf = data_status(newp, oldp)
                if ds == 'EOL':
                    eol = True
                    hash_map[hn] = ho
                elif ds == 'DIFF':
                    st = 'DIFFER'
            if j['log']:
                newp = os.path.join(outdir, 'new_logs', j['log'])
                ok, against, diff = (compare_logs(newp, [os.path.join(repo, 'results', x) for x in [j['log']] + j['alt_logs']],
                                                  KNOWN_DIFFS.get(j['name']), hash_map)
                                     if os.path.exists(newp) else (False, None, ['missing']))
                if ok == 'KNOWN':
                    if st == 'MATCH':
                        st = 'MATCH-KNOWN'
                elif not ok:
                    st = 'DIFFER'
                    open(os.path.join(outdir, 'new_logs', 'DIFF_' + j['log']), 'w', encoding='utf-8').write('\n'.join(diff) + '\n')
            if st == 'MATCH' and eol:
                st = 'MATCH-EOL'
        all_ok = all_ok and st in ('MATCH', 'MATCH-KNOWN', 'MATCH-EOL')
        lines.append(f'{st:<8} {rec["job"]}')
    lines.append('OVERALL (recompare): ' + ('ALL MATCH' if all_ok else 'NOT ALL MATCH'))
    print('\n'.join(lines))
    open(os.path.join(outdir, 'RECOMPARE.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    return 0 if all_ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--repo', required=True, help='folder rarita-schwinger-off-slice (from the phase1 zip)')
    ap.add_argument('--parent', help='folder rarita-schwinger-well-posedness at v1.0.11 (for step 1)')
    ap.add_argument('--out', default=None, help='results folder (default: results_offslice_<date_time>)')
    ap.add_argument('--skip-long', action='store_true', help='skip the ~12-minute dispersion proof')
    ap.add_argument('--timeout', type=int, default=5400, help='per-job timeout in seconds (default 5400)')
    ap.add_argument('--recompare', metavar='RESULTS_DIR',
                    help='do not run anything: re-compare the logs and data in an earlier results folder')
    a = ap.parse_args()
    if a.recompare:
        return recompare(os.path.abspath(a.repo), os.path.abspath(a.recompare))

    repo = os.path.abspath(a.repo)
    if not os.path.isfile(os.path.join(repo, 'SHA256SUMS.txt')) or not os.path.isdir(os.path.join(repo, 'scripts')):
        print(f'NOT VERIFIED: {repo} does not look like rarita-schwinger-off-slice (no SHA256SUMS.txt or scripts/).')
        return 2
    parent = os.path.abspath(a.parent) if a.parent else None
    if parent and not os.path.isdir(os.path.join(parent, 'scripts')):
        print(f'NOT VERIFIED: --parent {parent} does not look like rarita-schwinger-well-posedness (no scripts/).')
        return 2

    stamp = _dt.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    outdir = os.path.abspath(a.out or f'results_offslice_{stamp}')
    os.makedirs(outdir, exist_ok=True)
    work = os.path.join(outdir, 'work', 'rarita-schwinger-off-slice')
    shutil.copytree(repo, work)
    regen_dir = os.path.join(outdir, 'regenerated_data')
    os.makedirs(regen_dir, exist_ok=True)
    logs_dir = os.path.join(outdir, 'new_logs')
    os.makedirs(logs_dir, exist_ok=True)

    try:
        import sympy  # noqa
        sympy_v = sympy.__version__
    except Exception:
        sympy_v = 'NOT INSTALLED'
    try:
        import mpmath  # noqa
        mpmath_v = mpmath.__version__
    except Exception:
        mpmath_v = 'NOT INSTALLED'

    lines = []

    def say(s=''):
        print(s, flush=True)
        lines.append(s)

    say(f'{PKG_VERSION}')
    say(f'started {stamp} local | python {platform.python_version()} | sympy {sympy_v} | mpmath {mpmath_v} | {platform.system()} {platform.release()}')
    say(f'repo copy: work/rarita-schwinger-off-slice (your original folder is not modified)')

    # 1. checksums of the copy (paths in SHA256SUMS.txt start with rarita-schwinger-off-slice/)
    sums_path = os.path.join(work, 'SHA256SUMS.txt')
    say(f'SHA256SUMS.txt sha256 {sha256_file(sums_path)}')
    root = os.path.dirname(work)
    bad, n = [], 0
    for line in open(sums_path, encoding='utf-8'):
        line = line.strip()
        if not line:
            continue
        h, rel = line.split(None, 1)
        rel = rel.lstrip('*')
        p = os.path.join(root, rel.replace('/', os.sep))
        n += 1
        if not os.path.isfile(p) or sha256_file(p) != h:
            bad.append(rel)
    sums_ok = not bad
    say(f'checksums: {n - len(bad)}/{n} OK' + ('' if sums_ok else f'  MISMATCH: {bad[:10]}'))

    if sympy_v != '1.14.0':
        say(f'WARNING: SymPy {sympy_v}; the committed logs were made with 1.14.0. Differences may follow.')

    jobs = build_jobs(parent, a.skip_long)
    results = []
    committed_dir = os.path.join(work, 'results')
    t_all = time.time()
    for i, j in enumerate(jobs, 1):
        rec = {'job': j['name'], 'argv': None, 'status': None, 'committed_log': j['log']}
        if j.get('skip'):
            rec['status'] = 'NOT RUN'
            rec['reason'] = j['skip']
            say(f'[{i:2d}/{len(jobs)}] NOT RUN   {j["name"]}  ({j["skip"]})')
            results.append(rec)
            continue
        new_log = os.path.join(logs_dir, (j['log'] or (re.sub(r'\W+', '_', j['name']) + '.txt')))
        argv = []
        for x in j['argv']:
            if x == '@LOG':
                argv.append(new_log)
            elif x.startswith('@REGEN:'):
                argv.append(os.path.join(regen_dir, x[len('@REGEN:'):]))
            else:
                argv.append(x)
        cmd = [sys.executable, '-I'] + argv
        rec['argv'] = ['python', '-I'] + [os.path.relpath(x, work) if os.path.isabs(x) else x for x in argv]
        say(f'[{i:2d}/{len(jobs)}] running   {j["name"]} ...')
        t0 = time.time()
        stdout_path = os.path.join(logs_dir, 'stdout_' + re.sub(r'\W+', '_', j['name']) + '.txt')
        try:
            with open(stdout_path, 'w', encoding='utf-8') as so:
                cp = subprocess.run(cmd, cwd=work, stdout=so, stderr=subprocess.STDOUT, timeout=a.timeout)
            rc = cp.returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
        dt = time.time() - t0
        rec['seconds'] = round(dt, 1)
        rec['exit'] = rc
        if rc != 0:
            rec['status'] = 'FAIL'
            say(f'           FAIL      exit {rc} after {dt:.0f} s (see new_logs/{os.path.basename(stdout_path)})')
            results.append(rec)
            continue
        status = 'MATCH'
        notes = []
        hash_map, eol = {}, False
        for r in j['regen']:
            newp, oldp = os.path.join(regen_dir, r), os.path.join(work, 'data', r)
            if not os.path.exists(newp):
                status = 'DIFFER'
                notes.append(f'{r} not regenerated')
                continue
            ds, hn, ho, hf = data_status(newp, oldp)
            rec.setdefault('data', {})[r] = {'new': hn, 'committed': ho, 'result': ds, 'new_after_crlf_to_lf': hf}
            if ds == 'EOL':
                eol = True
                hash_map[hn] = ho
                notes.append(f'data/{r}: identical to the committed file except Windows line endings (CRLF); after CRLF->LF sha256 = committed')
            elif ds == 'DIFF':
                status = 'DIFFER'
                notes.append(f'data/{r} regenerated with a different SHA-256 (not only line endings)')
        if j['log']:
            if not os.path.exists(new_log):
                status = 'DIFFER'
                notes.append('script did not write its log')
            else:
                kd = KNOWN_DIFFS.get(j['name'])
                ok, against, diff = compare_logs(new_log, [os.path.join(committed_dir, j['log'])] +
                                                 [os.path.join(committed_dir, x) for x in j['alt_logs']], kd, hash_map)
                rec['log_sha256'] = sha256_file(new_log)
                rec['compared_against'] = against
                if ok == 'KNOWN':
                    if status == 'MATCH':
                        status = 'MATCH-KNOWN'
                    notes.append('only known, explained difference: ' + kd['why'])
                elif not ok:
                    status = 'DIFFER'
                    dpath = os.path.join(logs_dir, 'DIFF_' + j['log'])
                    open(dpath, 'w', encoding='utf-8').write('\n'.join(diff) + '\n')
                    notes.append(f'log differs from {against}: see new_logs/DIFF_{j["log"]}')
        if status == 'MATCH' and eol:
            status = 'MATCH-EOL'
        rec['status'] = status
        if notes:
            rec['notes'] = notes
        say(f'           {status:<9} {dt:.0f} s' + (('  ' + '; '.join(notes)) if notes else ''))
        results.append(rec)

    counts = {}
    for r in results:
        counts[r['status']] = counts.get(r['status'], 0) + 1
    all_match = sums_ok and all(r['status'] in ('MATCH', 'MATCH-KNOWN', 'MATCH-EOL') for r in results)
    say('')
    say('SUMMARY')
    say(f'  checksums: {"OK" if sums_ok else "MISMATCH"}')
    for k in ('MATCH', 'MATCH-KNOWN', 'MATCH-EOL', 'DIFFER', 'FAIL', 'NOT RUN'):
        if counts.get(k):
            say(f'  {k:<8} {counts[k]}')
    say(f'  total time {time.time() - t_all:.0f} s')
    if all_match:
        say('  OVERALL: ALL MATCH. The off-slice scripts reproduce their committed outputs on this machine'
            + (' (MATCH-KNOWN / MATCH-EOL lines are explained above).' if (counts.get('MATCH-KNOWN') or counts.get('MATCH-EOL')) else '.'))
    else:
        say('  OVERALL: NOT ALL MATCH. See the lines above; nothing has been changed. Send SUMMARY.txt and the DIFF_ files to Claude.')
    say('  Scope: reproducibility of software outputs only. Not a proof of correctness, not independent review.')

    manifest = {'package': PKG_VERSION, 'started_local': stamp, 'python': platform.python_version(),
                'sympy': sympy_v, 'mpmath': mpmath_v, 'os': f'{platform.system()} {platform.release()}',
                'sha256sums_ok': sums_ok, 'sha256sums_file_sha256': sha256_file(sums_path),
                'parent_given': bool(parent), 'skip_long': a.skip_long, 'jobs': results,
                'overall': 'ALL MATCH' if all_match else 'NOT ALL MATCH'}
    open(os.path.join(outdir, 'MANIFEST.json'), 'w', encoding='utf-8').write(json.dumps(manifest, indent=2) + '\n')
    open(os.path.join(outdir, 'SUMMARY.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
    # remove the work copy to keep the results folder small (inputs are re-derivable from the zip)
    shutil.rmtree(os.path.join(outdir, 'work'), ignore_errors=True)
    print(f'\nResults folder: {outdir}')
    return 0 if all_match else 1


if __name__ == '__main__':
    sys.exit(main())
