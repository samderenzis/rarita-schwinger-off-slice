"""m2_family_q_direction.py -- Milestone 2 (partial): constrained Jordan chains at BOTH off-slice crossings for EVERY
eps with 0 < |eps| <= 1, q = (1 + eps) q_slice, at the witness m, lam, b1, b2. Exact (SymPy), rational inequalities.

Usage:  python3 -I m2_family_q_direction.py [--out LOG] [--json data/m2_minor.json]

Argument (each item is checked below unless marked as imported):
  (R) Rotation invariance (m2_rotation_invariance.py): P(R^T xi) = Shat P(xi) Shat^-1 and L(R^T xi) Shat = S L(xi) for
      rotations in the (1,2) and (3,4) planes; with homogeneity, every xi with xi_b != 0 is similar to a point of the ray
      xi = (s, 0, 1, 0), s >= 0.  [imported]
  (M) Algebraic multiplicity (step 4, exact with eps symbolic): charpoly(P|ker L) along the ray is
      (x^2 - s^2 - 1)^4 Q(x^2, s^2, eps^2) up to a nonzero factor, so at a crossing (Q(1+s^2, s^2, eps^2) = 0) the
      eigenvalue r = |xi| has algebraic multiplicity >= 5 on ker L.  [imported]
  (C) Crossings: for 0 < eps^2 <= 1 the crossing polynomial C(u, eps^2) has two distinct positive roots u_-(eps) < u_+(eps)
      (discriminant and coefficient signs, below). They are continuous in eps and equal u* = 400/9801 at eps = 0.
  (G) Geometric multiplicity: with s = 2 tau/(1 - tau^2) (so |xi| = (1 + tau^2)/(1 - tau^2) is rational in tau) and
      eta = eps sqrt(2450894), the matrix M(tau, eta) = Dn * (P(xi) - |xi|) K(tau) (K: a basis of ker L(xi), Dn: a nonzero
      polynomial clearing denominators) has a fixed 8 x 8 minor d. With tau = 1/10 + a (tau = 1/10 <-> s* = 20/99), eta = b:
      on the box |a| <= ha, |b| <= hb,  sum_{(i,j) != 0} |c_ij| ha^i hb^j < |d(0,0)|  (exact rational arithmetic; |c| is
      bounded by |Re c| + |Im c|), so d != 0 and rank (P - r)|ker L >= rank M >= 8 there, i.e. geometric multiplicity <= 4.
  (B) Both crossings stay in the box for |eps| <= eps0: the crossing polynomial Ct(tau, eta) has no real root eta in
      [-hb, hb] on the box sides tau = 1/10 +- ha (exact Sturm count), and Ct(1/10, 0) = 0; by continuity of the two
      crossing branches tau_+-(eta) (which start at tau = 1/10) they cannot leave the box.
  Conclusion: geometric multiplicity <= 4 < 5 <= algebraic multiplicity at each crossing, so P|ker L is not
  diagonalisable at r: there is a nontrivial Jordan chain v -> w in ker L. The pi'-factorisation (step 2) makes this hold
  for every pi'. det S^0 != 0 for |eps| < 8.82 (step 2). Curved-symbol identification for all k and q:
  m2_curved_symbol_general.py. Coercivity: by reading the proof of Prop. 4.2 with q free.
Scope: a ONE-parameter family (q varies; m, lam, b1, b2 fixed at the witness values). For each eps the obstruction is a
separate pointwise statement with its own radius; no uniformity in eps is claimed. Not an open set in parameter space.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
from sympy import Rational as Q, sqrt  # noqa: E402
from sympy.polys.matrices import DomainMatrix  # noqa: E402
import rs_model as rm  # noqa: E402
import points as pts  # noqa: E402

ROWS = [2, 3, 4, 5, 8, 9, 10, 11]          # fixed choice (found by pivoting at the slice crossing; any nonzero minor works)
COLS = [0, 1, 4, 5, 8, 9, 10, 11]
EPS0 = [Q(1, 1000), Q(1, 100), Q(1, 10), Q(1, 2), Q(1)]
HA = {Q(1, 1000): Q(1, 20000), Q(1, 100): Q(1, 2000), Q(1, 10): Q(1, 200), Q(1, 2): Q(1, 50), Q(1): Q(1, 20)}
D = 2450894


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); ap.add_argument('--json'); args = ap.parse_args()
    log = []; ok_all = True; T0 = time.time()

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}  [{time.time() - T0:.1f} s]')

    say('Milestone 2 (partial) -- Jordan chains at both crossings for every 0 < |eps| <= 1 (witness m, lam, b; q free)')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__} | {platform.system()} {platform.machine()}')
    say('command: python3 -I scripts/m2_family_q_direction.py')
    for f in ['rs_model.py', 'points.py', 'm2_family_q_direction.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    tau, eta, a, b, v = sp.symbols('tau eta a b v', real=True)
    s = 2 * tau / (1 - tau ** 2); x = (1 + tau ** 2) / (1 - tau ** 2)

    say('C. Crossings are real, distinct and positive for 0 < eps^2 <= 1')
    c2, c1, c0 = [sp.expand(z) for z in pts.crossing_poly_coeffs(v)]
    disc = sp.factor(c1 ** 2 - 4 * c2 * c0)
    say(f'   discriminant in u, as a function of v = eps^2: {disc}')
    rec('discriminant > 0 on (0, 1]: v*(91809 - 823 v) > 0 there, with a positive constant factor',
        sp.expand(disc - sp.Integer(95227628033791228445188816896) * v * (91809 - 823 * v)) == 0)
    rec('C coefficient of u^2 > 0 on [0, 1] (no root there, positive at 0)', sp.Poly(c2, v).count_roots(0, 1) == 0 and c2.subs(v, 0) > 0)
    rec('C coefficient of u^0 > 0 on [0, 1]', sp.Poly(c0, v).count_roots(0, 1) == 0 and c0.subs(v, 0) > 0)
    rec('C coefficient of u^1 < 0 on [0, 1]  (so both roots u are positive)', sp.Poly(c1, v).count_roots(0, 1) == 0 and c1.subs(v, 0) < 0)

    say('G. The matrix M(tau, eta) and the minor d')
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam); W = rm.witness_model(q=qsl + eta / 35736)
    rec('tau = 1/10 gives s = 20/99 and |xi| = 101/99 (the slice crossing)',
        sp.simplify(s.subs(tau, Q(1, 10)) - Q(20, 99)) == 0 and sp.simplify(x.subs(tau, Q(1, 10)) - Q(101, 99)) == 0)
    P = W.Psym([s, 0, 1, 0]).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z))))
    ssym = sp.Symbol('s', real=True)
    Ks = sp.Matrix.hstack(*W.Lsym([ssym, 0, 1, 0]).nullspace()).applyfunc(sp.expand)     # basis of ker L, entries a + b*s
    rec('K(s) entries are affine in s', all(sp.Poly(z, ssym).degree() <= 1 for z in Ks if z != 0))
    L = W.Lsym([s, 0, 1, 0]); K = Ks.subs(ssym, s)
    rec('L K == 0 (columns of K lie in ker L)', (L * K).applyfunc(sp.cancel).is_zero_matrix)
    Dn = (tau ** 2 - 1) ** 2 * (4 * eta ** 2 - 27625 ** 2) * (16 * eta ** 2 - 63361 ** 2)
    M = ((P - x * sp.eye(16)) * K * Dn).applyfunc(sp.cancel)
    rec('Dn * (P - |xi|) K is polynomial in tau, eta', all(sp.denom(z) == 1 for z in M))
    sub = M.extract(ROWS, COLS).subs({tau: Q(1, 10) + a, eta: b}).applyfunc(sp.expand)
    d = DomainMatrix.from_Matrix(sub).det()
    dpoly = sp.Poly(sp.expand(DomainMatrix.from_Matrix(sub).domain.to_sympy(d)), a, b)
    terms = dpoly.terms(); d00 = dpoly.coeff_monomial(1)
    say(f'   minor rows {ROWS}, cols {COLS}: degree (a, b) = ({dpoly.degree(a)}, {dpoly.degree(b)}), {len(terms)} terms')
    say(f'   d(0,0) = {d00}')
    rec('d(0,0) is a nonzero rational', d00 != 0 and sp.im(d00) == 0)
    absb = lambda c: abs(sp.re(c)) + abs(sp.im(c))
    coef = [(m, absb(c)) for m, c in terms if m != (0, 0)]

    say('B + G. Box certificates')
    tst = Q(1, 10)
    u = sp.Symbol('u')
    Cuv = c2.subs(v, eta ** 2 / D) * u ** 2 + c1.subs(v, eta ** 2 / D) * u + c0.subs(v, eta ** 2 / D)
    Ct = sp.numer(sp.together(Cuv.subs(u, s ** 2)))
    Ct = sp.Poly(sp.expand(Ct), tau, eta).as_expr()
    rec('Ct(1/10, 0) == 0 (the centre is on the crossing curve)', Ct.subs({tau: tst, eta: 0}) == 0)
    certs = []
    for e0 in EPS0:
        hb = sp.ceiling(e0 * sqrt(D) * 1000) / Q(1000)
        ha = HA[e0]
        hb_ok = hb ** 2 >= e0 ** 2 * D
        inside = all(sp.Poly(Ct.subs(tau, tb), eta).count_roots(-hb, hb) == 0 for tb in (tst - ha, tst + ha))
        dn_ok = hb < Q(27625, 2) and tst + ha < 1 and tst - ha > -1        # Dn != 0 on the box
        bnd = sum(c * ha ** i * hb ** j for (i, j), c in coef)
        cert = bnd < abs(d00)
        ratio = bnd / abs(d00)
        certs.append({'eps0': str(e0), 'hb': str(hb), 'ha': str(ha), 'ratio': str(sp.N(ratio, 8))})
        rec(f'eps0 = {e0}: hb = {hb} (hb^2 >= eps0^2 D: {hb_ok}), ha = {ha}; crossings cannot leave the box: {inside}; '
            f'Dn != 0: {dn_ok}; Taylor bound / |d(0,0)| = {sp.N(ratio, 6)} < 1', hb_ok and inside and dn_ok and cert)
    say('   -> for every 0 < |eps| <= 1, rank (P - r)|ker L >= 8 at both crossings: geometric multiplicity <= 4 < 5.')

    if args.json:
        data = {'description': 'Milestone 2 minor: d(a, b) = det of rows/cols of Dn*(P - |xi|)K at tau = 1/10 + a, '
                               'eta = b; xi = (2tau/(1-tau^2), 0, 1, 0); q = q_slice + eta/35736; K entries are '
                               '[a_re, a_im, b_re, b_im] meaning (a + b*s), s = 2tau/(1-tau^2). Coefficients: [i, j, re, im].',
                'rows': ROWS, 'cols': COLS,
                'Dn': '(tau^2 - 1)^2 (4 eta^2 - 27625^2) (16 eta^2 - 63361^2)',
                'K': [[None] * 12 for _ in range(16)],
                'coefficients': [[m[0], m[1], str(sp.re(c)), str(sp.im(c))] for m, c in terms],
                'certificates': certs}
        for i in range(16):
            for j in range(12):
                pz = sp.Poly(Ks[i, j], ssym)
                c0_, c1_ = pz.coeff_monomial(1), pz.coeff_monomial(ssym)
                data['K'][i][j] = [str(sp.re(c0_)), str(sp.im(c0_)), str(sp.re(c1_)), str(sp.im(c1_))]
        with open(args.json, 'w', encoding='utf-8') as fh:
            json.dump(data, fh, indent=0)
        say(f'   wrote {os.path.relpath(args.json)}  sha256 {hashlib.sha256(open(args.json, "rb").read()).hexdigest()}')
    say()
    say('Scope: one-parameter family (q only); pointwise in eps; uses (R) rotation invariance, (M) step 4, step 2, the general')
    say('curved-symbol identification, and coercivity by proof reading. AI-checked only.')
    line = 'MILESTONE 2 (q-direction, |eps| <= 1): ALL PASS' if ok_all else 'MILESTONE 2 (q-direction): SOME CHECK FAILED'
    say(line)
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
