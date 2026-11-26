"""step2_reduction_q_free.py -- PLAN step 2: re-verify the reduction with the charge q FREE (an independent symbol).

Usage:  python3 -I step2_reduction_q_free.py [--out FILE]

Each check repeats a parent check (or a part of the paper's Prop. 5.1) with q replaced by an independent
symbol q, using rs_model. Exact arithmetic (sympy). Parent references are to rarita-schwinger-well-posedness
v1.0.11 and paper version 2.

  1. Clifford intertwining (parent check_1): involves m, lam, F only. q does not occur in X, Y, so it is
     q-independent by construction; this script confirms that X, Y contain no q.
  2. Divergence identity (paper Prop. 3.2, parent check_3): (D+Yhat)_A R^A[psi] = S^rho psi_rho on Model E,
     arbitrary psi, symbolic m, lam, b1, b2 AND q.
  3. Constraint intertwining (paper Prop. 3.4, parent check_4): Ktilde has no d/dt, and
     R^0[g] - star R^0[psi] + S^rho psi_rho = 0, arbitrary psi_0, Phi, symbolic m, lam, b1, b2, pi' AND q.
  4. pi'-factorisation (paper Prop. 5.1(a)): P_pi'(xi) - P_pi''(xi) = (pi' - pi'') M^{-1} T L(xi) for a
     symbolic covector xi, symbolic pi', pi'', symbolic q, at the witness m, lam, b1, b2 (paper M = MMmat).
  5. Invariance of ker L (paper Prop. 5.1(b)): L(xi) P(xi) = C_1(xi) L(xi), C_1 = pi' sum_a xi_a Gamma^0 Gamma^a,
     symbolic xi, symbolic pi', symbolic q, witness m, lam, b1, b2.
  6. q-independent pieces: L(xi) and M contain no q; L L^dag = 3|xi|^2 I for symbolic xi; det M != 0.
  7. det S^0 as a function of q at the witness m, lam, b1, b2: exact polynomial, its value on the slice,
     and its real roots (where the elimination of psi_0 fails). Roots are reported, not interpreted.

Checks 4 and 5 hold wherever det S^0 != 0 (S^0 is inverted); they are identities of rational functions in q
and xi, verified by expanding the difference to 0 after clearing nothing (sympy cancel/expand).
Scope: algebraic identities only. Coercivity and the curved-symbol identification are NOT checked here.
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
from sympy import Rational as Q, Matrix, zeros, eye  # noqa: E402
import rs_model as rm  # noqa: E402

qs = sp.Symbol('q', real=True)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=None); args = ap.parse_args()
    log = []

    def say(s=''):
        print(s, flush=True); log.append(s)

    ok_all = True
    say('PLAN step 2 -- reduction re-verified with q as an independent symbol')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__} | {platform.system()} {platform.machine()}')
    say('command: python3 -I scripts/step2_reduction_q_free.py')
    for f in ['rs_model.py', 'step2_reduction_q_free.py']:
        say(f'   {sha(os.path.join(HERE, f))}  scripts/{f}')
    say()

    def rec(name, good, t0):
        nonlocal ok_all
        ok_all &= bool(good)
        say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)} ({time.time() - t0:.1f} s)')

    # 1 ---------------------------------------------------------------------------------------------
    say('1. Clifford intertwining is q-independent')
    t0 = time.time(); Ms = rm.symbolic_model(q=qs)
    rec('X_A, Yhat_A contain no q (check_1 identity does not involve q)',
        not any(x.has(qs) for Mx in Ms.X + Ms.Y for x in Mx), t0)
    rec('Model.q is the free symbol', Ms.q == qs and not Ms.on_slice, t0)

    # 2 ---------------------------------------------------------------------------------------------
    say('2. Divergence identity, arbitrary psi, symbolic m, lam, b1, b2, q')
    t0 = time.time(); psi = rm.generic_psi(); R = Ms.R(psi)
    res = Ms.div(R) - sum((Ms.S[r] * psi[r] for r in rm.R5), zeros(4, 1))
    rec('S^rho contains q', any(x.has(qs) for Sx in Ms.S for x in Sx), t0)
    rec('no Float atoms', rm.no_float(*Ms.S, *Ms.X, *Ms.Y, *Ms.Om), t0)
    rec('DI holds (exact zero test on the chart)', all(rm.vanishes_on_chart(x) for x in res), t0)
    say(f'      zero-test statistics so far: {dict(rm.STATS)}')

    # 3 ---------------------------------------------------------------------------------------------
    say('3. Constraint intertwining, arbitrary psi_0, Phi, symbolic m, lam, b1, b2, pi\', q')
    t0 = time.time(); g = Ms.g_of(R)
    Kt = Matrix.vstack(*[psi[k].diff(rm.t) for k in rm.SP]) - g
    rec('Ktilde free of d/dt', all(not sp.expand(x).has(sp.Derivative(psi[r][c], rm.t))
                                   for x in Kt for r in rm.R5 for c in range(4)), t0)
    gpsi = [zeros(4, 1)] + [g[4 * (k - 1):4 * k, 0] for k in rm.SP]
    res = Ms.R(gpsi)[0] - Ms.star(R[0]) + sum((Ms.S[r] * psi[r] for r in rm.R5), zeros(4, 1))
    rec('CI holds (exact zero test on the chart)', all(rm.vanishes_on_chart(x) for x in res), t0)
    say(f'      zero-test statistics, checks 2+3: {dict(rm.STATS)}')
    say(f'      denominators cleared in step 2 of the zero test: {sorted(map(str, rm.DENOMS)) or "none"}')

    # symbolic covector and witness m, lam, b1, b2 with q free -----------------------------------------
    a1, a2, c1_, c2_ = sp.symbols('xi1 xi2 xi3 xi4', real=True); xs = [a1, a2, c1_, c2_]
    p1, p2 = sp.symbols('pip1 pip2')
    Wq = rm.witness_model(q=qs)
    T = rm.blk({(a - 1, 0): rm.G[0] * rm.G[a] for a in rm.SP}, 4, 1)          # (T rho)_a = Gamma^0 Gamma^a rho
    MMi = Wq.MMmat.inv()

    def P_at(pp):
        Wq.pip = pp
        return Wq.Psym(xs)

    # 4 ---------------------------------------------------------------------------------------------
    say('4. pi\'-factorisation P_pi\' - P_pi\'\' = (pi\' - pi\'\') M^-1 T L  (symbolic xi, pi\', pi\'\', q; witness m, lam, b)')
    t0 = time.time(); L = Wq.Lsym(xs)
    P1, P2 = P_at(p1), P_at(p2)
    diff = (P1 - P2 - (p1 - p2) * MMi * T * L).applyfunc(lambda e: sp.cancel(sp.expand(e)))
    rec('identity holds exactly', diff.is_zero_matrix, t0)
    rec('P depends on q (q is live in the symbol)', any(sp.expand(x).has(qs) for x in P1), t0)

    # 5 ---------------------------------------------------------------------------------------------
    say('5. ker L invariant: L P = C_1 L, C_1 = pi\' sum_a xi_a Gamma^0 Gamma^a  (symbolic xi, pi\', q; witness m, lam, b)')
    t0 = time.time()
    C1 = p1 * sum((xs[a - 1] * rm.G[0] * rm.G[a] for a in rm.SP), rm.Z)
    d5 = (L * P1 - C1 * L).applyfunc(lambda e: sp.cancel(sp.expand(e)))
    rec('identity holds exactly', d5.is_zero_matrix, t0)

    # 6 ---------------------------------------------------------------------------------------------
    say('6. q-independent pieces')
    t0 = time.time()
    rec('L(xi) contains no q', not any(x.has(qs) for x in L), t0)
    rec('M (= MMmat) contains no q', not any(x.has(qs) for x in Wq.MMmat), t0)
    rec('L L^dag == 3|xi|^2 I (symbolic real xi)',
        (L * L.H - 3 * sum(x * x for x in xs) * eye(4)).applyfunc(sp.expand).is_zero_matrix, t0)
    dM = Wq.MMmat.det(); rec(f'det M = {dM} != 0', dM != 0, t0)

    # 7 ---------------------------------------------------------------------------------------------
    say('7. det S^0 as a function of q (witness m, lam, b1, b2, pi\' irrelevant)')
    t0 = time.time()
    d = sp.factor(sp.expand(Wq.S[0].det()))
    say(f'      det S^0(q) = {d}')
    qslice = rm.pauli_q(Wq.m, Wq.lam)
    dval = sp.nsimplify(sp.radsimp(sp.expand(d.subs(qs, qslice))))
    rec(f'det S^0 at q = slice equals parent value 3063716808343140625/132101812606622109696 (got {dval})',
        dval == Q(3063716808343140625, 132101812606622109696), t0)
    num, den = sp.fraction(sp.together(sp.expand(d)))
    say(f'      denominator (independent of q): {sp.factor(den)}  | has q: {den.has(qs)}')
    numP = sp.Poly(sp.expand(num), qs)
    say(f'      numerator: degree {numP.degree()} in q; factor list: {sp.factor_list(sp.expand(num), qs)}')
    say(f'      q_slice = {sp.radsimp(qslice)} = {sp.N(qslice, 25)}')
    say('      complex roots of the numerator (exploration only, 30-digit numerics; ratio to q_slice):')
    for z in numP.nroots(n=30, maxsteps=200):
        say(f'         q = {sp.N(z, 25)}   q/q_slice = {sp.N(z / qslice, 20)}')
    say('7b. det S^0 in the parametrisation q = (1 + eps) q_slice (exact)')
    t0 = time.time(); eps, u = sp.symbols('epsilon u', real=True)
    de = sp.expand(sp.radsimp(sp.expand(d.subs(qs, (1 + eps) * qslice))))
    say(f'      det S^0(eps) coefficients (eps^4 .. eps^0): {sp.Poly(de, eps).all_coeffs()}')
    rec('det S^0(eps) has rational coefficients', all(c.is_Rational for c in sp.Poly(de, eps).all_coeffs()), t0)
    rec('det S^0(eps) is even in eps', sp.expand(de - de.subs(eps, -eps)) == 0, t0)
    fu = sp.factor(sp.expand(de.subs(eps, sp.sqrt(u))))
    say(f'      with u = eps^2: det S^0 = {fu}')
    us = sorted(sp.solve(fu, u))
    say(f'      zeros at eps^2 = {us} ~ {[sp.N(x, 15) for x in us]}')
    target = Q(27625**2, 9803576)
    rec(f'smallest zero is eps^2 = 27625^2/9803576 (|eps| ~ {sp.N(sp.sqrt(target), 15)}); det S^0 != 0 for all |eps| below it',
        us[0] == target and all(x > 0 for x in us), t0)
    say()
    say('Scope: algebraic reduction identities with q free. Not checked here: coercivity on a ball, the')
    say('curved-symbol identification, any crossing or Jordan structure off the slice.')
    line = 'STEP 2 (q free): ALL PASS' if ok_all else 'STEP 2 (q free): SOME CHECK FAILED'
    say(line)
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
