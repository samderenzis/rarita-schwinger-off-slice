"""m2_endpoint_gaps.py -- close the gaps 2.66 <= |eps| <= eps_L (lower crossing) and 6.28 <= |eps| <= eps_U (upper crossing),
endpoints included. Exact inequalities (fractions); SymPy only to build the polynomial matrices.

Usage:  python3 -I m2_endpoint_gaps.py [--out LOG] [--json data/m2_gap_certificate.json]

Lower gap: ray 1, xi = (s, 0, 1, 0), s = 2tau/(1 - tau^2), matrix M1 = Dn (P - |xi|) K (as in m2_family_extended.py).
  Ordinary segments, then one END segment with box [0, t_hi] in tau whose eps-range reaches past eps_L. The lower branch
  tau_-(eps) >= 0 is continuous up to eps_L, where it reaches tau = 0 (xi = (0, 0, 1, 0)).
Upper gap: ray 2, xi' = (1, 0, t, 0), t = 2sigma/(1 - sigma^2), |xi'| = (1 + sigma^2)/(1 - sigma^2). xi' is a positive
  multiple of (s, 0, 1, 0) with s = 1/t, so the Jordan structure at |xi'| equals that at |xi|. Basis K' of ker L(xi'):
  identity on the components outside block 2, block 2 solved through B(t) = Gamma^0 (Gamma^1 + t Gamma^3) Gamma^2,
  det B = (1 + t^2)^2 != 0, so K' is a basis for every real t, including t = 0. Matrix M2 = Dn' (P - |xi'|) K' with
  Dn' = (sigma^2 - 1)(sigma^2 + 1)^2 (4eta^2 - 27625^2)(16eta^2 - 63361^2). The upper crossing is the SMALLEST root
  sigma in (0, 1); its END segment has box [0, sigma_hi] and reaches past eps_U (xi' -> (1, 0, 0, 0)).
Each box: exact containment (m2_box_tools.segment_containment with ray and touch_zero) and the exact rank-8 bound
||C A - I||_inf < 1 (m2_box_tools.box_rank8). Both signs of eps.
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
import mpmath as mp  # noqa: E402
import rs_model as rm  # noqa: E402
import m2_box_tools as bt  # noqa: E402
from m2_family_extended import build_M, sqrt_lo, sqrt_hi, pivot, ratc  # noqa: E402

D = bt.D
VL = (sp.Integer(34404) - 12 * sp.sqrt(5659689)) / 823       # eps_L^2
VU = (sp.Integer(58228) - 20 * sp.sqrt(1645345)) / 823       # eps_U^2


def build_M2():
    sg, eta = sp.symbols('sigma eta', real=True); tsym = sp.Symbol('t', real=True)
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam); W = rm.witness_model(q=qsl + eta / 35736)
    Lt = W.Lsym([1, 0, tsym, 0]); Pc = [4, 5, 6, 7]; Fc = [i for i in range(16) if i not in Pc]
    B = Lt.extract(range(4), Pc)
    assert sp.factor(B.det()) == (tsym ** 2 + 1) ** 2
    Kt = sp.zeros(16, 12)
    for jj, f in enumerate(Fc):
        Kt[f, jj] = 1
    Bk = (-B.inv() * Lt.extract(range(4), Fc)).applyfunc(sp.cancel)
    for a_, p_ in enumerate(Pc):
        for jj in range(12):
            Kt[p_, jj] = Bk[a_, jj]
    assert (Lt * Kt).applyfunc(sp.cancel).is_zero_matrix
    t = 2 * sg / (1 - sg ** 2); x = (1 + sg ** 2) / (1 - sg ** 2)
    P = W.Psym([1, 0, t, 0]).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z))))
    K = Kt.subs(tsym, t).applyfunc(sp.cancel)
    Dn = (sg ** 2 - 1) * (sg ** 2 + 1) ** 2 * (4 * eta ** 2 - 27625 ** 2) * (16 * eta ** 2 - 63361 ** 2)
    M = ((P - x * sp.eye(16)) * K * Dn).applyfunc(sp.cancel)
    out = []
    for i in range(16):
        row = []
        for j in range(12):
            assert sp.denom(M[i, j]) == 1
            p = sp.Poly(sp.expand(M[i, j]), sg, eta)
            row.append([[a, b, str(sp.re(c)), str(sp.im(c))] for (a, b), c in p.terms()])
        out.append(row)
    # K'(t) * (1 + t^2): polynomial entries, recorded as coefficient lists [deg, re, im]
    Kj = [[[[k, str(sp.re(c)), str(sp.im(c))] for (k,), c in sp.Poly(sp.expand(sp.cancel(Kt[i, j] * (1 + tsym ** 2))), tsym).terms()]
           for j in range(12)] for i in range(16)]
    return out, Kj


def branch_param(e, branch, ray):
    """Numerical parameter (tau on ray 1, sigma on ray 2) of the branch at eps = e (planning only)."""
    c2, c1, c0 = [mp.mpf(z.numerator) / z.denominator for z in bt.c_coeffs_in_v(Fr(e) ** 2)]
    ds = c1 * c1 - 4 * c2 * c0
    u = [(-c1 - mp.sqrt(ds)) / (2 * c2), (-c1 + mp.sqrt(ds)) / (2 * c2)][0 if branch == 'lower' else 1]
    if u <= 0:
        return mp.mpf(0)
    s = mp.sqrt(u) if ray == 1 else 1 / mp.sqrt(u)
    return (mp.sqrt(1 + s * s) - 1) / s


def certify(Mpoly, Mnum, ray, branch, sign, e_lo, e_hi, touch_zero):
    ts = [branch_param(e_lo + (e_hi - e_lo) * n / 8, branch, ray) for n in range(9)]
    tmin, tmax = min(ts), max(ts); m = (tmax - tmin) * mp.mpf('0.25') + mp.mpf('1e-7')
    t_hi = Fr(int(mp.ceil((tmax + m) * 10 ** 10)), 10 ** 10)
    t_lo = Fr(0) if touch_zero else Fr(int(mp.floor((tmin - m) * 10 ** 10)), 10 ** 10)
    v_lo, v_hi = e_lo * e_lo * D, e_hi * e_hi * D
    if not bt.segment_containment(t_lo, t_hi, v_lo, v_hi, branch=branch, ray=ray, touch_zero=touch_zero):
        return None
    eta_a, eta_b = sqrt_lo(v_lo), sqrt_hi(v_hi)
    if sign < 0:
        eta_a, eta_b = -eta_b, -eta_a
    tc, ec = (t_lo + t_hi) / 2, (eta_a + eta_b) / 2; ha, hb = (t_hi - t_lo) / 2, (eta_b - eta_a) / 2
    A = mp.matrix(Mnum(mp.mpf(tc.numerator) / tc.denominator, mp.mpf(ec.numerator) / ec.denominator))
    rows, cols = pivot(A)
    Ci = mp.inverse(mp.matrix([[A[i, j] for j in cols] for i in rows]))
    C = [[ratc(Ci[i, j], den=2 ** 60 * 10 ** 20) for j in range(8)] for i in range(8)]
    ok, worst = bt.box_rank8(Mpoly, rows, cols, C, tc, ec, ha, hb)
    if not ok:
        return None
    return {'ray': ray, 'branch': branch, 'sign': sign, 'touch_zero': touch_zero, 'eps_lo': str(e_lo), 'eps_hi': str(e_hi),
            't_lo': str(t_lo), 't_hi': str(t_hi), 'v_lo': str(v_lo), 'v_hi': str(v_hi), 'eta_lo': str(eta_a), 'eta_hi': str(eta_b),
            'rows': rows, 'cols': cols, 'C': [[[str(c[0]), str(c[1])] for c in row] for row in C], 'row_sum_bound': float(worst)}


def cover(Mpoly, Mnum, ray, branch, sign, e_start, e_switch, e_end, say):
    boxes = []; e = e_start; step = Fr(1, 1000)
    while e < e_switch:
        h = min(step, e_switch - e)
        while True:
            bx = certify(Mpoly, Mnum, ray, branch, sign, e, e + h, False)
            if bx is not None:
                break
            h /= 2
            if h < Fr(1, 10 ** 9):
                raise SystemExit(f'cannot certify ray {ray} {branch} near eps = {float(e)}')
        boxes.append(bx); e += h; step = min(h * 2, Fr(1, 100))
    bx = certify(Mpoly, Mnum, ray, branch, sign, e_switch, e_end, True)
    if bx is None:
        raise SystemExit(f'end box failed: ray {ray} {branch} [{e_switch}, {e_end}]')
    boxes.append(bx)
    say(f'   ray {ray}, {branch} branch, eps {"+" if sign > 0 else "-"}: {len(boxes)} boxes covering [{e_start}, {e_end}] (end box touches 0)')
    return boxes


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); ap.add_argument('--json'); args = ap.parse_args()
    log = []; T0 = time.time(); ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    say('Endpoint gaps -- Jordan chains up to and including eps_L (lower) and eps_U (upper)')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__} | mpmath {mp.__version__}')
    say('command: python3 -I scripts/m2_endpoint_gaps.py')
    for f in ['rs_model.py', 'm2_box_tools.py', 'm2_family_extended.py', 'm2_endpoint_gaps.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    mp.mp.dps = 40
    eLs, eUs = Fr(26675, 10000), Fr(62913, 10000)          # rational endpoints of the end boxes
    rec(f'{eLs}^2 >= eps_L^2 and {eUs}^2 >= eps_U^2 (exact)',
        sp.Rational(eLs.numerator, eLs.denominator) ** 2 - VL >= 0 and sp.Rational(eUs.numerator, eUs.denominator) ** 2 - VU >= 0)
    say('Building M1 (ray 1) and M2 (ray 2) exactly ...')
    M1j, K1j = build_M(); M2j, K2j = build_M2()
    tau, eta = sp.symbols('tau eta')

    def num(Mj):
        Mp = [[bt.poly_from_json(e) for e in row] for row in Mj]
        Ms = sp.Matrix(16, 12, lambda i, j: sum((sp.Rational(c[0].numerator, c[0].denominator) + sp.I * sp.Rational(c[1].numerator, c[1].denominator)) * tau ** a * eta ** b
                                              for (a, b), c in Mp[i][j].items()))
        return Mp, sp.lambdify((tau, eta), Ms, 'mpmath')
    M1p, M1n = num(M1j); M2p, M2n = num(M2j)
    say(f'   done [{time.time() - T0:.1f} s]; max degrees of M2: sigma {max(i for r in M2p for p in r for i, j in p)}, eta {max(j for r in M2p for p in r for i, j in p)}')
    boxes = []
    for sign in (1, -1):
        boxes += cover(M1p, M1n, 1, 'lower', sign, Fr(266, 100), Fr(26674, 10000), eLs, say)
        boxes += cover(M2p, M2n, 2, 'upper', sign, Fr(628, 100), Fr(62911, 10000), eUs, say)
    worst = max(b['row_sum_bound'] for b in boxes)
    say(f'   total boxes: {len(boxes)}; largest row-sum bound {worst:.4f} < 1  [{time.time() - T0:.1f} s]')
    for ray, br, e0, e1 in [(1, 'lower', Fr(266, 100), eLs), (2, 'upper', Fr(628, 100), eUs)]:
        for sg_ in (1, -1):
            seg = [b for b in boxes if b['ray'] == ray and b['sign'] == sg_]
            cov = Fr(seg[0]['eps_lo']) == e0 and Fr(seg[-1]['eps_hi']) == e1 and seg[-1]['touch_zero'] and all(
                Fr(seg[i]['eps_hi']) == Fr(seg[i + 1]['eps_lo']) for i in range(len(seg) - 1))
            rec(f'ray {ray} {br}, sign {sg_:+d}: contiguous from {e0} to {e1}, last box touches 0', cov)
    if args.json:
        json.dump({'description': 'Endpoint-gap certificate. M1/K1: ray 1 (as m2_extended_certificate.json). M2: ray 2, terms '
                                  '[i, j, re, im] of sigma^i eta^j; K2: K\'(t)*(1+t^2) coefficient lists [deg, re, im].',
                   'M1': M1j, 'K1': K1j, 'M2': M2j, 'K2': K2j, 'boxes': boxes,
                   'end_eps': {'lower': str(eLs), 'upper': str(eUs)}}, open(args.json, 'w'))
        say(f'   wrote {os.path.relpath(args.json)}  sha256 {hashlib.sha256(open(args.json, "rb").read()).hexdigest()}')
    line = 'ENDPOINT GAPS: ALL PASS' if ok_all else 'ENDPOINT GAPS: SOME CHECK FAILED'
    say(line)
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
