"""coercivity_constants.py -- exact q-dependence of every constant in paper Prop. 4.2 (small-ball coercivity). SymPy.

Usage:  python3 -I audit/coercivity_constants.py [--out LOG]

Background as in the paper (b1 = 1, b2 = 1/3, chart V' with |theta_i - pi/2| < pi/3, |phi_i| < pi/3); m > 0, lam real
and q real are SYMBOLS here. R^0 Phi = sum_j L_j(y) d_j Phi + L_0(y) Phi is extracted from rs_model.Model.R (frame
formula of the paper), and M_0 = L_0^dag - sum_j mu^-1 d_j(mu L_j^dag), mu = sin th1 sin th2 / (b1^2 b2^2).
Checks / outputs:
  (1) R^0 contains no psi_0 and no time derivative; the decomposition is exact.
  (2) L_j contain no q (and no m, lam): the principal part of the constraint is charge-independent.
  (3) L L^dag = 3|xi|^2 I with L(xi) the frame symbol (symbolic xi): no q.
  (4) M_0 is affine in q, M_0 = M_0^(0) + q M_0^(1); exact Frobenius identity
        ||M_0||_F^2 = 3 c1^2 + c2^2/3 + 10 lam^2/3 + 16 m^2 + 12 q^2 (c1^2 + c2^2),  c_i = cot th_i,
      so on V' (|c_i| <= sqrt3):  m0(q) := sup ||M_0||_2 <= sup ||M_0||_F = sqrt(10 + 10 lam^2/3 + 16 m^2 + 72 q^2).
  (5) The other constants of Prop. 4.2 do not involve q: sigma = sqrt3 * min(b1, b2) = 1/sqrt3,
      l1 = 6 sqrt(b1^2 + b2^2) = 2 sqrt10, mu_- = 1/(4 b1^2 b2^2) = 9/4, Poincare constant 1/j_{1,1} < 1.
"""
import argparse, datetime, hashlib, os, platform, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..'); sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import sympy as sp  # noqa: E402
import rs_model as rm  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); args = ap.parse_args(); log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    say('Coercivity constants of paper Prop. 4.2 with q free (SymPy, exact)')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__}')
    for f in ['audit/coercivity_constants.py', 'scripts/rs_model.py']:
        say(f'   {hashlib.sha256(open(os.path.join(ROOT, f), "rb").read()).hexdigest()}  {f}')
    q = sp.Symbol('q', real=True); m, lam = sp.symbols('m lam', positive=True)
    b1, b2 = sp.Integer(1), sp.Rational(1, 3)
    W = rm.Model(m, lam, b1, b2, sp.Rational(1, 2), q=q)
    psi = rm.generic_psi(); R0 = W.R(psi)[0]
    th1, ph1, th2, ph2 = rm.th1, rm.ph1, rm.th2, rm.ph2; coords = [th1, ph1, th2, ph2]
    Phi = [psi[k][c] for k in range(1, 5) for c in range(4)]
    L0 = sp.zeros(4, 16); Lj = [sp.zeros(4, 16) for _ in range(4)]
    for i in range(4):
        ex = sp.expand(R0[i])
        for n, f in enumerate(Phi):
            for jj, y in enumerate(coords):
                Lj[jj][i, n] = ex.coeff(sp.Derivative(f, y))
            L0[i, n] = sp.expand(ex.subs({sp.Derivative(f, y): 0 for y in coords + [rm.t]})).coeff(f)
    reb = [sum(Lj[j][i, n] * sp.Derivative(Phi[n], coords[j]) for j in range(4) for n in range(16))
           + sum(L0[i, n] * Phi[n] for n in range(16)) for i in range(4)]
    rec('(1) R^0 = sum_j L_j d_j + L_0 exactly; no psi_0, no d/dt',
        all(sp.simplify(sp.expand(R0[i] - reb[i])) == 0 for i in range(4))
        and not any(sp.expand(R0[i]).has(psi[0][c]) or sp.expand(R0[i]).has(sp.Derivative(Phi[0], rm.t)) for i in range(4) for c in range(4)))
    rec('(2) L_j contain no q, m, lam', not any(z.has(q) or z.has(m) or z.has(lam) for M_ in Lj for z in M_))
    xs = sp.symbols('xi1:5', real=True)
    Lxi = W.Lsym(list(xs))
    rec('(3) L(xi) L(xi)^dag == 3|xi|^2 I and L(xi) contains no q',
        (Lxi * Lxi.H - 3 * sum(z * z for z in xs) * sp.eye(4)).applyfunc(sp.expand).is_zero_matrix and not any(z.has(q) for z in Lxi))
    mu = sp.sin(th1) * sp.sin(th2) / (b1 ** 2 * b2 ** 2)
    M0 = (L0.H - sum(((1 / mu) * sp.diff(mu * Lj[j].H, coords[j]) for j in range(4)), sp.zeros(16, 4))).applyfunc(sp.simplify)
    c1, c2 = sp.symbols('c1 c2', real=True)
    M0c = M0.applyfunc(lambda z: sp.expand(z.rewrite(sp.cot).subs({sp.cot(th1): c1, sp.cot(th2): c2})))
    rec('    M_0 depends on y only through c_i = cot(theta_i)', set().union(*[z.free_symbols for z in M0c]) <= {c1, c2, q, m, lam})
    rec('(4) M_0 is affine in q', all(sp.degree(z, q) <= 1 for z in M0c))
    F2 = sp.expand(sum(sp.re(z) ** 2 + sp.im(z) ** 2 for z in M0c))
    target = 3 * c1 ** 2 + c2 ** 2 / 3 + sp.Rational(10, 3) * lam ** 2 + 16 * m ** 2 + 12 * q ** 2 * (c1 ** 2 + c2 ** 2)
    rec('    ||M_0||_F^2 == 3c1^2 + c2^2/3 + 10 lam^2/3 + 16 m^2 + 12 q^2 (c1^2 + c2^2)', sp.expand(F2 - target) == 0)
    rec('    |cot theta| <= sqrt3 on [pi/6, 5pi/6] (cot decreasing; cot(pi/6) = sqrt3)', sp.cot(sp.pi / 6) == sp.sqrt(3))
    bound = sp.sqrt(10 + sp.Rational(10, 3) * lam ** 2 + 16 * m ** 2 + 72 * q ** 2)
    say(f'   => m0(q) <= {bound}')
    wit = {m: sp.sqrt(sp.Rational(823, 142944)), lam: sp.sqrt(sp.Rational(1, 2))}
    bw = sp.simplify(bound.subs(wit)); qsl = sp.sqrt(sp.Rational(823, 428832))
    say(f'   witness m, lam: m0(q) <= {bw};  on the slice q = {qsl}: {sp.N(bw.subs(q, qsl), 8)}')
    sigma = sp.sqrt(3) * sp.Min(b1, b2); l1 = 6 * sp.sqrt(b1 ** 2 + b2 ** 2); mum = 1 / (4 * b1 ** 2 * b2 ** 2)
    say(f'(5) q-free constants: sigma = {sigma} ~ {sp.N(sigma, 6)}, l1 = {sp.simplify(l1)} ~ {sp.N(l1, 6)}, mu_- = {mum}')
    r0 = sp.Min(sp.pi / 3, sigma / (2 * (l1 + bw)))
    for name, qq in [('slice', qsl), ('eps = 1', 2 * qsl), ('eps = 6.28', sp.Rational(728, 100) * qsl), ('q = 1', 1), ('q = 10', 10)]:
        say(f'   r0 >= {sp.N(r0.subs(q, qq), 6)} at {name}')
    c0 = (sigma / 2) * sp.sqrt(mum) / sp.sqrt(1 + (sp.pi / 3) ** 2)
    say(f'   c0 may be taken q-independent: c0 = (sigma/2) mu_-^(1/2) (1 + (pi/3)^2)^(-1/2) = {sp.N(c0, 6)} (valid for every r <= r0 <= pi/3)')
    say('(6) Appendix B constants (N1, N2 contain M_0 and its first derivatives)')
    rec('    M_0 is affine in c1 and c2 (no c1 c2, c_i^2 terms)', all(sp.Poly(z, c1, c2).total_degree() <= 1 for z in M0c if z != 0))
    A1 = M0c.applyfunc(lambda z: sp.diff(z, c1)); A2 = M0c.applyfunc(lambda z: sp.diff(z, c2))
    n1 = sp.expand(sum(sp.re(z) ** 2 + sp.im(z) ** 2 for z in A1)); n2 = sp.expand(sum(sp.re(z) ** 2 + sp.im(z) ** 2 for z in A2))
    rec(f'    ||dM_0/dc1||_F^2 == 3 + 12 q^2 (got {n1}) and ||dM_0/dc2||_F^2 == 1/3 + 12 q^2 (got {n2})',
        sp.expand(n1 - 3 - 12 * q ** 2) == 0 and sp.expand(n2 - sp.Rational(1, 3) - 12 * q ** 2) == 0)
    say('   d/d theta_i = -csc^2(theta_i) d/dc_i and csc^2 <= 4 on V\': ||d_theta1 M_0|| <= 4 sqrt(3 + 12 q^2), '
        '||d_theta2 M_0|| <= 4 sqrt(1/3 + 12 q^2), d_phi M_0 = 0.')
    say('   So N1 <= (q-free) + sqrt(10 + 10 lam^2/3 + 16 m^2 + 72 q^2), N2 <= (q-free) + 4 sqrt(3 + 12 q^2) + 4 sqrt(1/3 + 12 q^2):')
    say('   finite for every q, uniform on compact q-sets, growing linearly in |q| (not uniform on all of R).')
    line = 'COERCIVITY CONSTANTS: ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
