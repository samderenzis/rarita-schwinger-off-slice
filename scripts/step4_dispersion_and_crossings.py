"""step4_dispersion_and_crossings.py -- PLAN steps 3-4: q = (1 + eps) q_slice; where are the off-slice crossings?

Usage:  python3 -I step4_dispersion_and_crossings.py [--out LOG]

Witness m, lam, b1, b2, pi' = 1/2; covector ray xi = (s, 0, 1, 0). Parts marked EXACT are exact SymPy computations
(they are identities in s, eps); parts marked EXPLORATION are high-precision numerics (mpmath, 60 digits) and are
never results.

  A. EXACT  P depends on q only through eta = eps*sqrt(2450894) = 35736 (q - q_slice): P has entries in Q(i)(eta).
  B. EXACT  characteristic polynomial of P restricted to ker L(xi) (12 x 12), symbolic s and eps:
            (x^2 - s^2 - 1)^4 * Q(x^2, s^2, eps^2) / (det S^0 factors), Q quadratic in x^2, even in eps.
  C. EXACT  crossing polynomial C(u, eps) = Q(1 + u, u, eps^2), u = s^2; its discriminant in u, and in eps^2.
  D. EXPLORATION  for eps = 1/10, 1/100 (and 0): the two real crossings near s* = 20/99, and the singular values
            of (A - r)^j on ker L (A = P on ker L): numerically dim ker = 4, 5, 5 (one Jordan block) at each.
  E. EXACT  search: rational u on the conic w^2 = 4467 (u+1)(82875u+1267) near u* = 400/9801 with u, 1+u sums of two
            rational squares -> eps^2 rational, covector and eigenvalue rational. Reproduces points A and B (points.py).
"""
import argparse
import datetime
import hashlib
import os
import platform
import sys
import time
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
import mpmath as mp  # noqa: E402
from sympy import Rational as Q, sqrt  # noqa: E402
from sympy.ntheory import factorint  # noqa: E402
import rs_model as rm  # noqa: E402
import points as pts  # noqa: E402


def norm_qi(x):
    """Positive rational that is a sum of two rational squares."""
    x = Q(x)
    return x > 0 and all(not (p % 4 == 3 and e % 2) for p, e in list(factorint(x.p).items()) + list(factorint(x.q).items()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); args = ap.parse_args()
    log = []; T0 = time.time()

    def say(s=''):
        print(s, flush=True); log.append(s)
    say('PLAN steps 3-4 -- off-slice dispersion relation and crossings near the witness')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__} | mpmath {mp.__version__}')
    say('command: python3 -I scripts/step4_dispersion_and_crossings.py')
    for f in ['rs_model.py', 'points.py', 'step4_dispersion_and_crossings.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    s, e, u, v, eta = sp.symbols('s epsilon u v eta', real=True)
    x = sp.Symbol('x')   # plain symbol: charpoly(x).as_expr() returns an assumption-free x
    D = 2450894
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam)

    say('A. EXACT: dependence on q')
    We = rm.witness_model(q=qsl + eta / 35736)
    Pe = We.Psym([s, 0, 1, 0]).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z))))
    rad = any(a.is_Pow and a.exp.is_Rational and not a.exp.is_integer for z in Pe for a in sp.preorder_traversal(z))
    say(f'   P(s,0,1,0) with q = q_slice + eta/35736 has no radical (entries in Q(i)(s, eta)): {not rad}')

    say('B. EXACT: characteristic polynomial of P on ker L(xi), xi = (s, 0, 1, 0), symbolic eps')
    W = rm.witness_model(q=(1 + e) * qsl); xi = [s, 0, 1, 0]
    P = W.Psym(xi).applyfunc(lambda z: sp.radsimp(sp.expand(z)))
    L = W.Lsym(xi); Kb = sp.Matrix.hstack(*L.nullspace())
    A = ((Kb.H * Kb).applyfunc(sp.simplify).inv() * Kb.H * P * Kb).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z))))
    inv_ok = (P * Kb - Kb * A).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z)))).is_zero_matrix
    say(f'   P ker L subset ker L (P K = K A exactly): {inv_ok}')
    cp = sp.cancel(A.charpoly(x).as_expr())
    n, d = sp.fraction(sp.together(cp)); n = sp.expand(sp.radsimp(sp.expand(n)))
    say(f'   no radical, no i in the numerator: {not n.has(sqrt(D)) and not n.has(sp.I)}')
    fl = sp.factor_list(n, x, s, e)[1]
    light = [f for f, mlt in fl if mlt == 4]; Qp = [f for f, mlt in fl if mlt == 1][0]
    say(f'   light-branch factor (x^2 - s^2 - 1) with multiplicity 4: {len(light) == 1 and sp.expand(light[0]**2 - (x**2 - s**2 - 1)**2) == 0}')
    say(f'   Q(x, s, eps) = {Qp}')
    say(f'   denominator: {sp.factor(d)}')
    say(f'   Q even in eps: {sp.expand(Qp - Qp.subs(e, -e)) == 0}')
    say(f'   Q at eps = 0: {sp.factor(Qp.subs(e, 0))}   (the parent double G-sheet)')
    say(f'   ({time.time() - T0:.1f} s)')

    say('C. EXACT: crossings x^2 = |xi|^2 = 1 + s^2 of the light branch with the Q-sheets')
    C = sp.expand(Qp.subs(x, sqrt(1 + u)).subs(s, sqrt(u)))
    say(f'   C(u, eps) coefficients in u: {[sp.factor(c) for c in sp.Poly(C, u).all_coeffs()]}')
    say(f'   C(u, 0) = {sp.factor(C.subs(e, 0))}   (double root u* = 400/9801, s* = 20/99)')
    disc_u = sp.factor(sp.discriminant(C, u))
    say(f'   discriminant in u: {disc_u}')
    say('   -> two distinct real roots u for 0 < eps^2 < 91809/823 (~ 10.56^2); split ~ linear in eps')
    Cv = sp.Poly(sp.expand(C.subs(e, sqrt(v))), v); a2, a1, a0 = Cv.all_coeffs()
    disc_v = sp.factor(a1**2 - 4 * a2 * a0)
    say(f'   as a quadratic in v = eps^2: discriminant {disc_v}')
    say(f'   ({time.time() - T0:.1f} s)')

    say('D. EXPLORATION (60-digit numerics; not a result): Jordan structure at the crossings')
    mp.mp.dps = 60

    def svals(M):
        return sorted(abs(z) for z in mp.svd_c(M, compute_uv=False))
    for ev in [Q(1, 10), Q(1, 100), Q(0)]:
        Wn = rm.witness_model(q=(1 + ev) * qsl)
        roots = [Q(400, 9801)] if ev == 0 else sp.Poly(C.subs(e, ev), u).nroots(n=60)
        for ur in roots:
            sv = Q(20, 99) if ev == 0 else sp.sqrt(sp.Float(ur, 60))
            Pn, Ln = Wn.Psym([sv, 0, 1, 0]), Wn.Lsym([sv, 0, 1, 0])
            Pm = mp.matrix([[mp.mpmathify(sp.N(Pn[i, j], 60)) for j in range(16)] for i in range(16)])
            Lm = mp.matrix([[mp.mpmathify(sp.N(Ln[i, j], 60)) for j in range(16)] for i in range(4)])
            _, _, Vh = mp.svd_c(Lm, full_matrices=True)
            Kn = mp.matrix(16, 12)
            for c in range(12):
                for i in range(16):
                    Kn[i, c] = mp.conj(Vh[4 + c, i])
            An = Kn.H * Pm * Kn; r = mp.sqrt(1 + mp.mpf(sp.N(ur, 60)))
            Nn = An - r * mp.eye(12); out = []
            Mj = mp.eye(12)
            for j in range(1, 4):
                Mj = Mj * Nn; sv_ = svals(Mj); out.append(sum(1 for z in sv_ if z < mp.mpf(10) ** -40))
            say(f'   eps = {ev}: s = {mp.nstr(mp.sqrt(mp.mpf(sp.N(ur, 60))), 18)}, r = {mp.nstr(r, 18)}: '
                f'numerical dim ker (A - r)^j, j = 1..3: {out}  (|PK - KA| = {mp.nstr(mp.mnorm(Pm * Kn - Kn * An, 1), 2)})')
    say(f'   ({time.time() - T0:.1f} s)')

    say('E. EXACT search for a point with rational eps^2, covector and eigenvalue (reproduces point A)')
    F = 4467 * (u + 1) * (82875 * u + 1267)
    u0 = Q(400, 9801); w0 = sp.sqrt(F.subs(u, u0)); m_ = sp.Symbol('m')
    say(f'   conic w^2 = 4467 (u+1)(82875u+1267) through (u*, {w0}); 4467 = square-free part of 2^10 3^3 823^2 1489^5')
    poly = sp.Poly(sp.expand((w0 + m_ * (u - u0))**2 - F), u)
    c2, c1, c0 = poly.all_coeffs(); usec = sp.cancel(c0 / (c2 * u0))
    mt = sp.diff(F, u).subs(u, u0) / (2 * w0)
    say(f'   tangent slope at u*: {mt}; second intersections for slopes m = {mt} + k/den, den < 10, |k| <= 60')
    cands = []
    for den in range(1, 10):
        for k in range(-60, 61):
            if gcd(abs(int(mt) * den + k), den) != 1:
                continue
            uu = usec.subs(m_, Q(int(mt) * den + k, den))
            if uu == u0 or not (uu > 0 and abs(uu - u0) < Q(5, 1000)) or not (norm_qi(uu) and norm_qi(1 + uu)):
                continue
            A_, B_, C_ = a2.subs(u, uu), a1.subs(u, uu), a0.subs(u, uu)
            sq = sp.sqrt(B_**2 - 4 * A_ * C_)
            vv = min(z for z in [(-B_ + sq) / (2 * A_), (-B_ - sq) / (2 * A_)] if z > 0)
            cands.append((len(str(uu)) + len(str(vv)), uu, vv))
    cands.sort()
    for c in cands[:5]:
        say(f'   u = {c[1]}, eps^2 = {c[2]} (eps ~ {sp.N(sp.sqrt(c[2]), 8)})')
    for nm in ['A', 'B']:
        P_ = pts.POINTS[nm]; UU, VV = Q(P_['U'].numerator, P_['U'].denominator), Q(P_['V'].numerator, P_['V'].denominator)
        found = [i for i, c in enumerate(cands) if c[1] == UU and c[2] == VV]
        say(f'   point {nm} (points.py; {P_["BRANCH"]} crossing) is candidate number {found[0] + 1 if found else None} by size; '
            f'U {"<" if UU < u0 else ">"} u*')
    say(f'   ({time.time() - T0:.1f} s)')
    say('Scope: B, C, E are exact; D is exploration only. The points are checked independently in step 5.')
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
