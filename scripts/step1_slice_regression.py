"""step1_slice_regression.py -- PLAN step 1: the q-independent code must reproduce the parent exactly ON the slice.

Usage:
    python3 -I step1_slice_regression.py --parent <path to rarita-schwinger-well-posedness v1.0.11> [--out FILE]

What it does (all exact, sympy):
  0. Checks the SHA-256 of every parent file it uses against the v1.0.11 values pinned below.
  A. Runs the parent's five check scripts UNMODIFIED, each in a fresh process, with rs_model standing in
     for bhmap_defs and q passed explicitly as lam*sqrt(2/3)*m (parent_shim_runner.py). Compares stdout
     with the parent's golden_outputs.txt, section by section, character for character.
  B. Runs witness_check.py (the witness_min.py checks rebuilt on rs_model) and compares with the golden
     witness_min.py section (the final 'sympy x.y.z' line is omitted, as in the golden file).
  C. Matrix-level comparison at the witness: S^rho, P(xi*), L(xi*), M from rs_model (q explicit) against
     the parent's own witness_min.py and bhmap_defs.py objects, entry by entry (difference expands to 0).
  D. q really is free: with a symbol q, S^rho is affine in q, and substituting q -> lam*sqrt(2/3)*m gives
     the parent's symbolic S^rho exactly; at the witness parameters, a q off the slice changes S^rho.
Exit code 0 only if every comparison passes.

Scope: this is a regression test of the new code ON the slice (code consistency with the parent).
It is not an independent confirmation of anything, and it says nothing about points off the slice.
"""
import argparse
import contextlib
import datetime
import hashlib
import io
import os
import platform
import runpy
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

PINNED = {  # rarita-schwinger-well-posedness v1.0.11, from its SHA256SUMS
    'results/golden_outputs.txt': '290340280bbad09005580c7092cf97229576249e174d51e78939c7212e2ad851',
    'scripts/bhmap_defs.py': '3d0fcd3c3a0ed88aef4ec33e1a014d1621e744621541ca0b2ded5ed6dd446bc0',
    'scripts/check_1_clifford_intertwining.py': '2a627df9a3b407c03ce9c4926b1c6d9d7e3ffce3a8623f9166c14606b42a4ddb',
    'scripts/check_2_geometry.py': '81470aca217fb276823a59d9778ff17274f6b8cdbd69c2083bafc206fcac0568',
    'scripts/check_3_divergence_identity.py': '1cbd75bcadaa8d6d398ca922355a3aa6c2c1cee67bc0bbad28f2c5ec3b75f836',
    'scripts/check_4_constraint_intertwining.py': '722c80750f13af2173216a23ef0f8d6d72494408a138c08e85b0f7144d5c703b',
    'scripts/check_5_principal_symbol.py': 'e1f20c0cce59a79ab4e093cf2d574bfd2e71068f86ce7e894b94fab9dc3cf0e6',
    'scripts/witness_min.py': '4c79837ac4d04fd0255ff5b999930c162f69d37b4d2933972ea29a0f7f87bc64',
}
CHECKS = ['check_1_clifford_intertwining.py', 'check_2_geometry.py', 'check_3_divergence_identity.py',
          'check_4_constraint_intertwining.py', 'check_5_principal_symbol.py']
OURS = ['rs_model.py', 'parent_shim_runner.py', 'witness_check.py', 'step1_slice_regression.py']


def sha(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def golden_sections(path):
    """Parse golden_outputs.txt into {script name: [lines]} (blank lines dropped)."""
    secs, cur = {}, None
    for line in open(path, encoding='utf-8').read().splitlines():
        if line.startswith('=== '):
            cur = line[4:].strip(); secs[cur] = []
        elif cur is not None and line.strip():
            secs[cur].append(line)
    return secs


def run(cmd):
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--parent', required=True)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    parent = os.path.abspath(args.parent)
    log = []

    def say(s=''):
        print(s); log.append(s)

    ok_all = True
    import sympy as sp
    say('PLAN step 1 -- slice regression of rs_model (q independent) against the parent, v1.0.11')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__} | {platform.system()} {platform.machine()}')
    say(f'command: python3 -I scripts/step1_slice_regression.py --parent <parent repo>')
    say()

    say('0. Parent input files (SHA-256 against pinned v1.0.11 values)')
    for rel, h in PINNED.items():
        got = sha(os.path.join(parent, rel)); good = got == h; ok_all &= good
        say(f'   {"OK  " if good else "FAIL"} {got}  {rel}')
    say('   This repository\'s scripts used:')
    for f in OURS:
        say(f'        {sha(os.path.join(HERE, f))}  scripts/{f}')
    if not ok_all:
        say('Parent files do not match v1.0.11; stopping.'); return finish(log, args.out, False)
    gold = golden_sections(os.path.join(parent, 'results/golden_outputs.txt'))
    say()

    say('A. Parent check scripts, unmodified, run on rs_model with q passed explicitly')
    for c in CHECKS:
        rc, out, err, dt = run([sys.executable, '-I', os.path.join(HERE, 'parent_shim_runner.py'),
                                os.path.join(parent, 'scripts', c)])
        lines = [x for x in out.splitlines() if x.strip()]
        good = rc == 0 and lines == gold[c]; ok_all &= good
        say(f'   {"OK  " if good else "FAIL"} {c}: exit {rc}, {len(lines)} lines, identical to golden: {lines == gold[c]} ({dt:.1f} s)')
        if not good:
            say('      --- got ---'); [say('      ' + x) for x in lines]
            say('      --- stderr ---'); [say('      ' + x) for x in err.splitlines()[-15:]]
    say()

    say('B. witness_check.py (witness_min.py checks rebuilt on rs_model, q explicit)')
    rc, out, err, dt = run([sys.executable, '-I', os.path.join(HERE, 'witness_check.py')])
    lines = [x for x in out.splitlines() if x.strip() and not x.startswith('sympy ')]
    good = rc == 0 and lines == gold['witness_min.py']; ok_all &= good
    say(f'   {"OK  " if good else "FAIL"} witness_check.py: exit {rc}, {len(lines)} lines, identical to golden witness_min.py: '
        f'{lines == gold["witness_min.py"]} ({dt:.1f} s)')
    for x in lines:
        say('      ' + x)
    if not good:
        say('      --- stderr ---'); [say('      ' + x) for x in err.splitlines()[-15:]]
    say()

    say('C. Matrix-level comparison at the witness (entrywise, exact)')
    sys.path.insert(0, HERE)
    import rs_model
    with contextlib.redirect_stdout(io.StringIO()):
        pw = runpy.run_path(os.path.join(parent, 'scripts/witness_min.py'))
    pdefs_ns = runpy.run_path(os.path.join(parent, 'scripts/bhmap_defs.py'))
    m, lam = sp.sqrt(sp.Rational(823, 142944)), sp.sqrt(sp.Rational(1, 2))
    Wn = rs_model.Model(m, lam, 1, sp.Rational(1, 3), sp.Rational(1, 2), q=lam * sp.sqrt(sp.Rational(2, 3)) * m)
    Wp = pdefs_ns['witness_model']()
    xs = [sp.Rational(20, 99), 0, 1, 0]

    def same(A, B):
        return A.shape == B.shape and all(sp.expand(x) == 0 for x in (A - B))
    comps = [
        ('S^rho (5 matrices) vs parent witness_min.py', all(same(Wn.S[i], pw['S'][i]) for i in range(5))),
        ('S^rho (5 matrices) vs parent bhmap_defs.witness_model()', all(same(Wn.S[i], Wp.S[i]) for i in range(5))),
        ('X_A (5 matrices) vs parent witness_min.py', all(same(Wn.X[i], pw['X'][i]) for i in range(5))),
        ('paper M = MMmat vs parent witness_min.py MM', same(Wn.MMmat, pw['MM'])),
        ('Mmat = -(paper M) (sign convention kept)', same(Wn.Mmat, -pw['MM'])),
        ('P(xi*) vs parent witness_min.py Psym(xi*)', same(Wn.Psym(xs), pw['Psym'](xs))),
        ('P(xi*) vs parent bhmap_defs Psym(xi*)', same(Wn.Psym(xs), Wp.Psym(xs))),
        ('L(xi*) vs parent witness_min.py Lsym(xi*)', same(Wn.Lsym(xs), pw['Lsym'](xs))),
        ('q value equals parent witness_min.py q', sp.expand(Wn.q - pw['q']) == 0),
    ]
    for name, good in comps:
        ok_all &= good
        say(f'   {"OK  " if good else "FAIL"} {name}: {good}')
    say()

    say('D. q is a genuinely free parameter')
    qs = sp.Symbol('q')
    Ms = rs_model.symbolic_model(q=qs); Mp = pdefs_ns['symbolic_model']()
    affine = all(sp.degree(sp.expand(x), qs) <= 1 for Sx in Ms.S for x in Sx)
    has_q = any(sp.expand(x).has(qs) for Sx in Ms.S for x in Sx)
    back = all(same(Ms.S[i].subs(qs, Mp.lam * sp.sqrt(sp.Rational(2, 3)) * Mp.m), Mp.S[i]) for i in range(5))
    Wq = rs_model.witness_model(q=sp.Rational(11, 10) * Wn.q)
    moved = not all(same(Wq.S[i], Wn.S[i]) for i in range(5))
    dflt = rs_model.witness_model()
    dflt_ok = dflt.on_slice and all(same(dflt.S[i], Wn.S[i]) for i in range(5))
    for name, good in [('symbolic S^rho depends on q', has_q),
                       ('symbolic S^rho is affine (degree <= 1) in q', affine),
                       ('S^rho(q -> lam*sqrt(2/3)*m) == parent symbolic S^rho', back),
                       ('witness S^rho changes when q = (11/10) x slice value', moved),
                       ('default q=None gives the slice (same S^rho as explicit q)', dflt_ok)]:
        ok_all &= good
        say(f'   {"OK  " if good else "FAIL"} {name}: {good}')
    say()
    say('Scope: regression of the new code ON the slice against the parent code (code consistency only).')
    say('Nothing here is evaluated off the slice except the single check that S^rho changes with q.')
    return finish(log, args.out, ok_all)


def finish(log, out, ok):
    line = 'STEP 1 SLICE REGRESSION: ALL PASS' if ok else 'STEP 1 SLICE REGRESSION: FAILED'
    print(line); log.append(line)
    if out:
        with open(out, 'w', encoding='utf-8') as fh:
            fh.write('\n'.join(log) + '\n')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
