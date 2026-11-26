"""step5_curved_symbol.py -- curved-symbol identification at one named off-slice point (analogue of parent check_5).

Usage:  python3 -I step5_curved_symbol.py --point A [--out LOG]

At y0 = (pi/2, 0, pi/2, 0), the principal symbol of the CURVED eliminated operator K' (psi_0 = W Phi) in the
chart covector k equals i * P(xi) with xi = E(y0)^T k = (k1, k2, k3/3, k4/3) (b1 = 1, b2 = 1/3, sin th_i = 1).
Here xi = XI of the point, so k = (XI_1, XI_2, 3 XI_3, 3 XI_4). Plane wave Phi = c exp(i Lam k.y), c in C^16 symbolic,
exactly as in the parent's check_5. The charge is q = q_slice + eta/35736 with eta a SYMBOL, so the identity is
checked for every q at once (in particular at the point's eta = ETA_RAT sqrt(DELTA)), with the
witness m, lam, b1, b2 and pi' = 1/2. Exact; SymPy.
"""
import argparse
import datetime
import hashlib
import os
import platform
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
from sympy import I, Rational as Q, Matrix  # noqa: E402
import rs_model as rm  # noqa: E402
import points as pts  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--point', required=True, choices=list(pts.POINTS))
    ap.add_argument('--out'); args = ap.parse_args(); name = args.point; p = pts.POINTS[name]
    log = []; T0 = time.time()

    def say(s=''):
        print(s, flush=True); log.append(s)
    say(f'PLAN step 5 -- curved-symbol identification at point {name} (q symbolic through eta)')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__}')
    say(f'command: python3 -I scripts/step5_curved_symbol.py --point {name}')
    for f in ['rs_model.py', 'points.py', 'step5_curved_symbol.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    eta = sp.Symbol('eta')
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam)
    W = rm.witness_model(q=qsl + eta / 35736)
    xi = [Q(x.numerator, x.denominator) for x in p['XI']]
    kk = [xi[0] / W.b1, xi[1] / W.b1, xi[2] / W.b2, xi[3] / W.b2]          # sin th_i = 1 at y0
    Lam = sp.Symbol('Lam', real=True)
    S0i = W.S[0].inv(); Wm = [-S0i * W.S[k] for k in rm.SP]
    cs = sp.symbols('c0:16')
    ph = sp.exp(I * Lam * (kk[0] * rm.th1 + kk[1] * rm.ph1 + kk[2] * rm.th2 + kk[3] * rm.ph2))
    Phi = [Matrix(cs[4 * j:4 * j + 4]) * ph for j in range(4)]
    psi = [sum((Wm[k - 1] * Phi[k - 1] for k in rm.SP), sp.zeros(4, 1))] + Phi
    Kp = (-W.g_of(W.R(psi)) / ph).applyfunc(sp.expand)
    y0 = {rm.th1: sp.pi / 2, rm.ph1: 0, rm.th2: sp.pi / 2, rm.ph2: 0}
    lead = Kp.applyfunc(lambda z: z.coeff(Lam, 1).subs(y0))
    P = W.Psym(xi)
    deg_ok = all(sp.degree(z, Lam) <= 1 for z in Kp)
    ident = all(sp.cancel(sp.radsimp(sp.expand(z))) == 0 for z in (lead - I * P * Matrix(cs)))
    nofl = rm.no_float(lead, P)
    say(f'   k = E(y0)^-T XI = {tuple(str(z) for z in kk)}')
    say(f'   {"OK  " if nofl else "FAIL"} no Float atoms: {nofl}')
    say(f'   {"OK  " if deg_ok else "FAIL"} degree in Lam <= 1: {deg_ok}')
    say(f'   {"OK  " if ident else "FAIL"} principal symbol of curved K\' at y0 == i P(XI), for every q (eta symbolic): {ident}')
    say(f'   ({time.time() - T0:.1f} s)')
    ok = nofl and deg_ok and ident
    say(f'CURVED-SYMBOL IDENTIFICATION AT POINT {name}: ' + ('PASS' if ok else 'FAIL'))
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
