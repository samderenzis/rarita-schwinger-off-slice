"""step5_point_exact.py -- PLAN step 5 at one named off-slice point (see points.py). Exact; SymPy 1.14.0.

Usage:  python3 -I step5_point_exact.py --point A [--out LOG] [--json data/point_A_chain.json]

All linear algebra is done over the number field K = Q(i, sqrt(DELTA)) of the point, with SymPy's
DomainMatrix over QQ.algebraic_field(I, sqrt(DELTA)) (exact). Steps:
  0. points.self_check(name) (fractions).
  1. Build S^rho, P(XI), L(XI) from rs_model with q = q_slice + ETA/35736, ETA a symbol; confirm that no radical
     survives (entries in Q(i)(ETA)); then put ETA = ETA_RAT*sqrt(DELTA).
  2. det S^0 at A: nonzero, rational, equal to the step 2 formula det S^0(eps) at eps^2 = V.
  3. L L^dag = 3|XI|^2 I.
  4. Characteristic polynomial of P(XI) on C^16: rational coefficients, exact factorisation, multiplicity of R,
     reality of every root.
  5. On ker L(XI): dim ker (P - R)^j for j = 1..4.
  6. Chain: v spans (P - R) ker (P - R)^2 restricted to ker L (one-dimensional here), w solves (P - R) w = v in ker L
     (free variables set to 0); scaled to have integer coordinates over the basis 1, i, s, i*s.
     Checks: L v = 0, L w = 0, (P - R) v = 0, (P - R) w = v, v != 0; and the gamma-identity of the parent.
  7. pi'-independence: the chain relations hold for pi' in {0, 1/2, 2, -3/7} (paper Prop. 5.1(a), re-checked
     with q free in step 2).
  8. Writes the exact data (parameters, XI, R, v, w) to JSON for the independent fractions-only verifier.
Scope: algebra at ONE point (pointwise). Curved-symbol identification: step5_curved_symbol.py.
Coercivity: see the result note. Nothing here is a family statement.
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
from math import lcm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
from sympy import Rational as Q, sqrt, I  # noqa: E402
from sympy.polys.matrices import DomainMatrix  # noqa: E402
import rs_model as rm  # noqa: E402
import points as pts  # noqa: E402


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def Qr(f):
    return Q(f.numerator, f.denominator)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--point', required=True, choices=list(pts.POINTS))
    ap.add_argument('--out'); ap.add_argument('--json'); args = ap.parse_args()
    name = args.point; p = pts.POINTS[name]
    log = []; ok_all = True; T0 = time.time()

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good)
        say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}  [{time.time() - T0:.1f} s]')

    say(f'PLAN step 5 -- exact constrained Jordan structure at off-slice point {name}')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__} | {platform.system()} {platform.machine()}')
    say(f'command: python3 -I scripts/step5_point_exact.py --point {name}')
    for f in ['rs_model.py', 'points.py', 'step5_point_exact.py']:
        say(f'   {sha(os.path.join(HERE, f))}  scripts/{f}')
    say()

    say(f'0. Point {name} constants (fractions)')
    for k, v in pts.self_check(name).items():
        rec(k, v)
    s = sqrt(p['DELTA']); xi = [Qr(x) for x in p['XI']]; r = Qr(p['R'])
    say(f'   eps = {"+" if p["SIGN"] > 0 else "-"}sqrt({p["V"]}) ~ {sp.N(p["SIGN"] * sqrt(Qr(p["V"])), 15)};  q/q_slice = 1 + eps;  {p["BRANCH"]} crossing')
    say(f'   ETA = eps*sqrt(2450894) = {Qr(p['ETA_RAT'])}*sqrt({p['DELTA']});  number field K = Q(i, sqrt({p['DELTA']}))')
    say(f'   XI = {tuple(str(x) for x in xi)}')
    say(f'   R = |XI| = {r}')
    say()

    say('1. Symbols at A')
    eta = sp.Symbol('eta')
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam)
    W = rm.witness_model(q=qsl + eta / 35736)
    rec('q_slice == sqrt(2450894)/35736', sp.expand(qsl - sqrt(2450894) / 35736) == 0)
    S = [Sx.applyfunc(lambda z: sp.radsimp(sp.expand(z))) for Sx in W.S]
    P_eta = W.Psym(xi).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z))))
    no_rad = not any(a.is_Pow and a.exp.is_Rational and not a.exp.is_integer
                     for M in S + [P_eta] for z in M for a in sp.preorder_traversal(z))
    rec('S^rho and P(XI) have entries in Q(i)(ETA): no radical left', no_rad)
    eta_val = Qr(p['ETA_RAT']) * s
    Pk = P_eta.subs(eta, eta_val); L = W.Lsym(xi)
    K = sp.QQ.algebraic_field(I, s)
    say(f'   field: {K}')
    dP = DomainMatrix.from_Matrix(Pk).convert_to(K); dL = DomainMatrix.from_Matrix(L).convert_to(K)
    rec('L(XI) has entries in Q(i) and no q', not any(z.has(eta) for z in L))

    say('2. det S^0 at A')
    dS0 = DomainMatrix.from_Matrix(S[0].subs(eta, eta_val)).convert_to(K).det()
    dS0s = sp.nsimplify(sp.expand(K.to_sympy(dS0)))
    say(f'   det S^0 = {dS0s}')
    e2 = Qr(p['V'])
    formula = (9803576 * e2 - 763140625) * (39214304 * e2 - 4014616321) / Q(132101812606622109696)
    rec('det S^0 != 0', dS0s != 0)
    rec('det S^0 equals the step 2 formula at eps^2 = V', sp.expand(dS0s - formula) == 0)

    say('3. Constraint symbol')
    rec('L L^dag == 3|XI|^2 I', (L * L.H - 3 * sum(x * x for x in xi) * sp.eye(4)).applyfunc(sp.expand).is_zero_matrix)
    rec('rank L(XI) == 4 (dim ker L = 12)', dL.rank() == 4)

    say('4. Characteristic polynomial of P(XI) (16 x 16)')
    x = sp.Symbol('x')
    cp = [K.to_sympy(c) for c in dP.charpoly()]
    rational = all(sp.expand(c).is_Rational for c in cp)
    rec('coefficients rational', rational)
    cpx = sp.Poly([sp.expand(c) for c in cp], x)
    fl = sp.factor_list(cpx.as_expr())
    say(f'   factorisation: {fl}')
    mult_r = sum(m for f, m in fl[1] if sp.Poly(f, x).degree() == 1 and sp.solve(f, x)[0] == r)
    rec('R is a root of multiplicity 5', mult_r == 5)
    allreal = all(sp.Poly(f, x).degree() == 1 or
                  (sp.Poly(f, x).degree() == 2 and sp.discriminant(f, x) > 0) for f, m in fl[1])
    rec('all 16 roots real (linear factors, quadratic factors with positive discriminant)', allreal)

    say('5. Constrained structure on ker L(XI)')
    N = dP - DomainMatrix.eye(16, K) * K.convert(r)
    kb = dL.nullspace().transpose()
    rec('ker L basis: 16 x 12', kb.shape == (16, 12))
    rec('ker L invariant under P: rank [kb | P kb] == 12', kb.hstack(dP * kb).rank() == 12)
    dims = []; Nj = kb
    for j in range(1, 5):
        Nj = N * Nj; dims.append(12 - Nj.rank())
    say(f'   dim ker (P - R)^j on ker L, j = 1..4: {dims}')
    rec('ranks stabilised (dims[2] == dims[3])', dims[2] == dims[3])
    say(f'   algebraic multiplicity of R on ker L = {dims[-1]}; geometric multiplicity = {dims[0]}')
    positive = dims[1] > dims[0]
    outcome = 'POSITIVE (nontrivial Jordan chain in ker L)' if positive else 'NEGATIVE (R semisimple on ker L)'
    say(f'   OUTCOME at this crossing: {outcome}')
    if not positive:
        say('   no chain to extract')

    vi = wi = None; den = 1
    if positive:
        say('6. Explicit chain v -> w in ker L')
        Nk = N * kb
        Z2 = (N * Nk).nullspace()
        yv = None
        for i in range(Z2.shape[0]):
            y = Z2[i:i + 1, :].transpose()
            if not (Nk * y).is_zero_matrix:
                yv = y; break
        v = Nk * yv
        vl = v.to_list(); idx = next(i for i in range(16) if vl[i][0] != K.zero)
        inv = K.quo(K.one, vl[idx][0])
        v = DomainMatrix([[K.mul(vl[i][0], inv)] for i in range(16)], (16, 1), K)          # v[idx] = 1
        # particular w with free variables 0: rref of [Nk | v]
        Rm, piv = Nk.hstack(v).to_field().rref()
        assert 12 not in piv, 'N w = v not solvable in ker L'
        Rl = Rm.to_list(); yl = [[K.zero] for _ in range(12)]
        for row, pc in enumerate(piv):
            yl[pc][0] = Rl[row][12]
        w = kb * DomainMatrix(yl, (12, 1), K)
        sD = sqrt(p['DELTA'])

        def coords(e):
            e = sp.expand(e)
            a, b = e.coeff(sD, 0).coeff(I, 0), e.coeff(sD, 0).coeff(I, 1)
            c, d = e.coeff(sD, 1).coeff(I, 0), e.coeff(sD, 1).coeff(I, 1)
            assert all(z.is_Rational for z in (a, b, c, d)) and sp.expand(a + b * I + (c + d * I) * sD - e) == 0
            return [a, b, c, d]
        vc = [coords(K.to_sympy(z[0])) for z in v.to_list()]; wc = [coords(K.to_sympy(z[0])) for z in w.to_list()]
        den = lcm(*[int(c.q) for e in vc + wc for c in e])
        vi = [[int(c * den) for c in e] for e in vc]; wi = [[int(c * den) for c in e] for e in wc]
        from math import gcd
        g0 = 0
        for e in vi + wi:
            for c in e:
                g0 = gcd(g0, c)
        vi = [[c // g0 for c in e] for e in vi]; wi = [[c // g0 for c in e] for e in wi]
        scale = Q(den, g0)
        v = v * K.convert(scale); w = w * K.convert(scale)
        den = 1
        say(f'   v, w scaled by {scale} to coprime integer coordinates over 1, i, s, i*s (s = sqrt({p['DELTA']}))')
        say(f'   largest coordinate digits: v {max(len(str(abs(c))) for e in vi for c in e)}, w {max(len(str(abs(c))) for e in wi for c in e)}')
        rec('L v == 0', (dL * v).is_zero_matrix)
        rec('L w == 0', (dL * w).is_zero_matrix)
        rec('(P - R) v == 0', (N * v).is_zero_matrix)
        rec('(P - R) w == v', (N * w - v).is_zero_matrix)
        rec('v != 0', not v.is_zero_matrix)
        g = sp.Symbol('gamma', positive=True)
        vS, wS, NS = v.to_Matrix(), w.to_Matrix(), N.to_Matrix()
        aG = vS - I * g * wS
        rec('L a_gamma == 0 and (gamma - i(P-R)) a_gamma == -i gamma^2 w',
            (L * aG).applyfunc(sp.expand).is_zero_matrix and
            ((g * sp.eye(16) - I * NS) * aG + I * g ** 2 * wS).applyfunc(sp.expand).is_zero_matrix)

        say('7. Independence of pi\' (chain relations at other pi\')')
        for pp in [Q(0), Q(1, 2), Q(2), Q(-3, 7)]:
            W.pip = pp
            Pp = W.Psym(xi).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z)))).subs(eta, eta_val)
            Np = DomainMatrix.from_Matrix(Pp).convert_to(K) - DomainMatrix.eye(16, K) * K.convert(r)
            rec(f'pi\' = {pp}: (P-R)v == 0 and (P-R)w == v', (Np * v).is_zero_matrix and (Np * w - v).is_zero_matrix)
        W.pip = Q(1, 2)

    if args.json:
        data = {
            'point': name, 'branch': p['BRANCH'], 'sign_of_eps': p['SIGN'], 'eps_squared': str(p['V']),
            'ETA_RAT': str(p['ETA_RAT']),
            'description': f'Off-slice point {name}: exact constrained Jordan structure (pointwise). Entries of v, w are '
                           '[a, b, c, d] meaning (a + b*i + (c + d*i)*sqrt(DELTA)) / DEN.',
            'parameters': {'lam^2': '1/2', 'm^2': '823/142944', 'b1': '1', 'b2': '1/3', 'pi_prime': '1/2',
                           'q': 'q_slice*(1+eps), q_slice = sqrt(2450894)/35736, eps = SIGN*sqrt(eps_squared)',
                           'q_minus_q_slice': f"({p['ETA_RAT']})*sqrt({p['DELTA']})/35736"},
            'DELTA': str(p['DELTA']), 'XI': [str(x) for x in xi], 'R': str(r), 'DEN': str(den),
            'outcome': outcome,
            'v': None if vi is None else [[str(c) for c in e] for e in vi],
            'w': None if wi is None else [[str(c) for c in e] for e in wi],
            'ker_dims_P_minus_R_on_kerL': dims, 'det_S0': str(dS0s),
            'charpoly_factorisation': str(fl),
        }
        with open(args.json, 'w', encoding='utf-8') as fh:
            json.dump(data, fh, indent=1)
        say(f'   wrote {os.path.relpath(args.json)}  sha256 {sha(args.json)}')
    say()
    say('Scope: one parameter point, one covector, one eigenvalue (pointwise). Not a family statement.')
    line = f'STEP 5 POINT {name} (exact): ALL PASS' if ok_all else f'STEP 5 POINT {name} (exact): SOME CHECK FAILED'
    say(line)
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
