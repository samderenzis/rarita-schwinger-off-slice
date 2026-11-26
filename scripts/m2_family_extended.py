"""m2_family_extended.py -- extended Milestone 2 certificate: Jordan chains along the WHOLE real crossing set, as far as
certified, for both signs of eps. Exact inequalities (fractions); SymPy only to build the polynomial matrix.

Usage:  python3 -I m2_family_extended.py [--out LOG] [--json data/m2_extended_certificate.json]

Setting as m2_family_q_direction.py (witness m, lam, b1, b2; q = (1 + eps) q_slice; ray xi = (s, 0, 1, 0),
s = 2 tau/(1 - tau^2), eta = eps sqrt(2450894)). Exact facts used:
  * lower crossing u_-(eps) > 0 exactly for 0 < eps^2 < eps_L^2 = (34404 - 12 sqrt(5659689))/823 (eps_L ~ 2.6674),
    where the constant coefficient c0 of C(u, eps^2) changes sign (the direction degenerates to s = 0);
  * upper crossing u_+(eps) finite and positive exactly for 0 < eps^2 < eps_U^2 = (58228 - 20 sqrt(1645345))/823
    (eps_U ~ 6.2912), where the leading coefficient c2 vanishes (s -> infinity);
  * discriminant > 0 for 0 < eps^2 < 91809/823; det S^0 != 0 for |eps| < 8.82.
Certificate: |eps| <= 1 for both branches is m2_family_q_direction.py. Here each branch and each sign of eps is
covered from |eps| = 1 up to a target (lower 2.66, upper 6.28) by eps-segments; each segment has a box
[t_lo, t_hi] x [eta range] with (i) exact containment of that branch (m2_box_tools.segment_containment) and
(ii) an exact rank-8 certificate ||C A - I||_inf < 1 (m2_box_tools.box_rank8). Boxes are found by bisection.
Conclusion on each certified segment: geometric multiplicity <= 4 < 5 <= algebraic, so a Jordan chain in ker L.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import sys
import time
from fractions import Fraction as Fr
from math import isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
import mpmath as mp  # noqa: E402
import rs_model as rm  # noqa: E402
import m2_box_tools as bt  # noqa: E402

TARGET = {'lower': Fr(266, 100), 'upper': Fr(628, 100)}
D = bt.D


def build_M():
    tau, eta = sp.symbols('tau eta', real=True); ssym = sp.Symbol('s', real=True)
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam); W = rm.witness_model(q=qsl + eta / 35736)
    s = 2 * tau / (1 - tau ** 2); x = (1 + tau ** 2) / (1 - tau ** 2)
    P = W.Psym([s, 0, 1, 0]).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z))))
    Ks = sp.Matrix.hstack(*W.Lsym([ssym, 0, 1, 0]).nullspace())
    assert (W.Lsym([ssym, 0, 1, 0]) * Ks).applyfunc(sp.expand).is_zero_matrix
    K = Ks.subs(ssym, s)
    Dn = (tau ** 2 - 1) ** 2 * (4 * eta ** 2 - 27625 ** 2) * (16 * eta ** 2 - 63361 ** 2)
    M = ((P - x * sp.eye(16)) * K * Dn).applyfunc(sp.cancel)
    out = []
    for i in range(16):
        row = []
        for j in range(12):
            assert sp.denom(M[i, j]) == 1
            p = sp.Poly(sp.expand(M[i, j]), tau, eta)
            row.append([[a, b, str(sp.re(c)), str(sp.im(c))] for (a, b), c in p.terms()])
        out.append(row)
    Kj = []
    for i in range(16):
        r_ = []
        for j in range(12):
            pz = sp.Poly(sp.expand(Ks[i, j]), ssym); a0, a1 = pz.coeff_monomial(1), pz.coeff_monomial(ssym)
            r_.append([str(sp.re(a0)), str(sp.im(a0)), str(sp.re(a1)), str(sp.im(a1))])
        Kj.append(r_)
    return out, Kj


def sqrt_lo(v, den=10 ** 6):
    """Rational r <= sqrt(v) (v >= 0 Fraction)."""
    r = Fr(isqrt(int(v * den * den)), den)
    while r * r > v:
        r -= Fr(1, den)
    return r


def sqrt_hi(v, den=10 ** 6):
    r = sqrt_lo(v, den)
    while r * r < v:
        r += Fr(1, den)
    return r


def branch_tau(e, k):
    """Numerical tau of branch k (0 lower, 1 upper) at eps = e (planning only)."""
    c2, c1, c0 = [mp.mpf(z.numerator) / z.denominator for z in bt.c_coeffs_in_v(Fr(e) ** 2)]
    ds = c1 * c1 - 4 * c2 * c0
    u = [(-c1 - mp.sqrt(ds)) / (2 * c2), (-c1 + mp.sqrt(ds)) / (2 * c2)][k]
    s = mp.sqrt(u)
    return (mp.sqrt(1 + s * s) - 1) / s


def pivot(A):
    A = A.copy(); rows = list(range(16)); cols = list(range(12)); pr = []; pc = []
    for _ in range(8):
        best = (-1, None, None)
        for i in rows:
            for j in cols:
                if abs(A[i, j]) > best[0]:
                    best = (abs(A[i, j]), i, j)
        _, i0, j0 = best; pr.append(i0); pc.append(j0); rows.remove(i0); cols.remove(j0)
        for i in rows:
            f = A[i, j0] / A[i0, j0]
            for j in range(12):
                A[i, j] -= f * A[i0, j]
    return sorted(pr), sorted(pc)


def ratc(z, den=2 ** 60):
    return (Fr(round(float(mp.re(z)) * den), den) if abs(mp.re(z)) < 1e300 else Fr(0),
            Fr(round(float(mp.im(z)) * den), den))


def try_segment(Mpoly, Mnum, k, sign, e_lo, e_hi):
    """Return a certified box dict, or None."""
    ts = [branch_tau(e_lo + (e_hi - e_lo) * n / 8, k) for n in range(9)]
    tmin, tmax = min(ts), max(ts); m = (tmax - tmin) * mp.mpf('0.25') + mp.mpf('1e-6')
    t_lo = Fr(int(mp.floor((tmin - m) * 10 ** 9)), 10 ** 9); t_hi = Fr(int(mp.ceil((tmax + m) * 10 ** 9)), 10 ** 9)
    v_lo, v_hi = e_lo * e_lo * D, e_hi * e_hi * D
    if not bt.segment_containment(t_lo, t_hi, v_lo, v_hi, branch=['lower', 'upper'][k]):
        return None
    eta_a, eta_b = sqrt_lo(v_lo), sqrt_hi(v_hi)               # |eta| range, rational, enclosing
    if sign < 0:
        eta_a, eta_b = -eta_b, -eta_a
    tc, ec = (t_lo + t_hi) / 2, (eta_a + eta_b) / 2; ha, hb = (t_hi - t_lo) / 2, (eta_b - eta_a) / 2
    A = mp.matrix(Mnum(mp.mpf(tc.numerator) / tc.denominator, mp.mpf(ec.numerator) / ec.denominator))
    rows, cols = pivot(A)
    Asub = mp.matrix([[A[i, j] for j in cols] for i in rows])
    Ci = mp.inverse(Asub)
    C = [[ratc(Ci[i, j], den=2 ** 60 * 10 ** 20) for j in range(8)] for i in range(8)]
    ok, worst = bt.box_rank8(Mpoly, rows, cols, C, tc, ec, ha, hb)
    if not ok:
        return None
    return {'branch': ['lower', 'upper'][k], 'sign': sign, 'eps_lo': str(e_lo), 'eps_hi': str(e_hi),
            't_lo': str(t_lo), 't_hi': str(t_hi), 'v_lo': str(v_lo), 'v_hi': str(v_hi),
            'eta_lo': str(eta_a), 'eta_hi': str(eta_b), 'rows': rows, 'cols': cols,
            'C': [[[str(c[0]), str(c[1])] for c in row] for row in C], 'row_sum_bound': float(worst)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); ap.add_argument('--json'); args = ap.parse_args()
    log = []; T0 = time.time(); ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)
    say('Milestone 2 extended -- Jordan chains along the real crossing set (both branches, both signs of eps)')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__} | mpmath {mp.__version__}')
    say('command: python3 -I scripts/m2_family_extended.py')
    for f in ['rs_model.py', 'm2_box_tools.py', 'm2_family_extended.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    v = sp.Symbol('v')
    c2 = 677329 * v ** 2 - 95843288 * v + 2732361984; c0 = 677329 * v ** 2 - 56628984 * v + 368640000
    eL2 = (34404 - 12 * sp.sqrt(5659689)) / 823; eU2 = (58228 - 20 * sp.sqrt(1645345)) / 823
    say('Exact thresholds')
    for nm, poly, root in [('c0 (lower branch hits s = 0)', c0, eL2), ('c2 (upper branch reaches s = infinity)', c2, eU2)]:
        good = sp.simplify(poly.subs(v, root)) == 0 and min(sp.solve(poly, v)) == root
        ok_all &= good
        say(f'   {"OK  " if good else "FAIL"} smallest positive root of {nm}: eps^2 = {root} ~ {sp.N(root, 12)} (eps ~ {sp.N(sp.sqrt(root), 10)}): {good}')
    tl, tu = TARGET['lower'], TARGET['upper']
    good = (tl ** 2 < eL2) == True and (tu ** 2 < eU2) == True
    ok_all &= good
    say(f'   {"OK  " if good else "FAIL"} targets {tl} < eps_L and {tu} < eps_U: {good}')
    say('Building M(tau, eta) exactly ...')
    Mjson, Kjson = build_M(); Mpoly = [[bt.poly_from_json(e) for e in row] for row in Mjson]
    tau, eta = sp.symbols('tau eta')
    Msym = sp.Matrix(16, 12, lambda i, j: sum((sp.Rational(c[0].numerator, c[0].denominator) + sp.I * sp.Rational(c[1].numerator, c[1].denominator)) * tau ** a * eta ** b
                                            for (a, b), c in Mpoly[i][j].items()))
    Mnum = sp.lambdify((tau, eta), Msym, 'mpmath'); mp.mp.dps = 40
    say(f'   done [{time.time() - T0:.1f} s]')
    boxes = []
    for k, name in [(0, 'lower'), (1, 'upper')]:
        for sign in (1, -1):
            e = Fr(1); target = TARGET[name]; step = Fr(1, 10); nseg = 0
            while e < target:
                h = min(step, target - e)
                while True:
                    bx = try_segment(Mpoly, Mnum, k, sign, e, e + h)
                    if bx is not None:
                        break
                    h /= 2
                    if h < Fr(1, 10 ** 5):
                        raise SystemExit(f'cannot certify {name} sign {sign} near eps = {float(e)}')
                boxes.append(bx); nseg += 1; e += h
                step = min(h * 2, Fr(1, 2))
            say(f'   {name} branch, eps {"+" if sign > 0 else "-"}: segments covering [1, {target}]: {nseg}  [{time.time() - T0:.1f} s]')
    worst = max(b['row_sum_bound'] for b in boxes)
    say(f'   total boxes: {len(boxes)}; largest row-sum bound ||C A - I||_inf <= {worst:.4f} < 1')
    # coverage check (exact): consecutive, starting at 1, ending at target
    for name in TARGET:
        for sign in (1, -1):
            seg = [b for b in boxes if b['branch'] == name and b['sign'] == sign]
            cov = Fr(seg[0]['eps_lo']) == 1 and Fr(seg[-1]['eps_hi']) == TARGET[name] and all(
                Fr(seg[i]['eps_hi']) == Fr(seg[i + 1]['eps_lo']) for i in range(len(seg) - 1))
            ok_all &= cov
            say(f'   {"OK  " if cov else "FAIL"} {name}, sign {sign:+d}: segments contiguous from 1 to {TARGET[name]}: {cov}')
    say('   -> with m2_family_q_direction.py (|eps| <= 1): Jordan chains at the lower crossing for 0 < |eps| <= 2.66 and at')
    say('      the upper crossing for 0 < |eps| <= 6.28 (both signs).')
    if args.json:
        json.dump({'description': 'Extended Milestone 2 certificate. M entries: [i, j, re, im] terms of tau^i eta^j. '
                                  'Boxes: see scripts/m2_box_tools.py. K: [a_re, a_im, b_re, b_im] = a + b*s.', 'M': Mjson,
                   'K': Kjson, 'targets': {k: str(v) for k, v in TARGET.items()},
                   'boxes': boxes}, open(args.json, 'w', encoding='utf-8'))
        say(f'   wrote {os.path.relpath(args.json)}  sha256 {hashlib.sha256(open(args.json, "rb").read()).hexdigest()}')
    line = 'MILESTONE 2 EXTENDED: ALL PASS' if ok_all else 'MILESTONE 2 EXTENDED: SOME CHECK FAILED'
    say(line)
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
