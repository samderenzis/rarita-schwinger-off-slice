"""endpoint_multiplicity.py -- exact algebraic multiplicity at the two endpoint crossings (SymPy, exact).

Usage:  python3 -I audit/endpoint_multiplicity.py [--out LOG]

Q(x, s, eps) from data/m2_dispersion.json (proved: det(x - P(s,0,1,0)) = (x^2 - s^2 - 1)^4 (x^2 - (1+s^2)/4)^2 Q/lc_x(Q)).
 (1) Q has total degree <= 4 in (x, s); so Q_h(x, s, b) := b^4 Q(x/b, s/b, eps) is a polynomial, and by homogeneity of
     det(x - P(xi)) (degree 16 in (x, xi)) we get, for t != 0,
       det(x - P(1,0,t,0)) = (x^2 - 1 - t^2)^4 (x^2 - (1+t^2)/4)^2 Q_h(x, 1, t)/lc_x(Q_h),
     and by continuity in t also at t = 0. lc_x(Q_h) = lc_x(Q) (the x^4 coefficient does not involve s).
 (2) Lower endpoint: xi = (0, 0, 1, 0), eps^2 = eps_L^2 = (34404 - 12 sqrt(5659689))/823, r = 1:
     Q(1, 0, eps_L) = 0 and d/dX Q(X = 1) != 0 (X = x^2): r = 1 is a simple root of Q, so alg = 4 + 1 = 5.
 (3) Upper endpoint: xi = (1, 0, 0, 0), eps^2 = eps_U^2 = (58228 - 20 sqrt(1645345))/823, r = 1:
     Q_h(1, 1, 0) = 0 and d/dX Q_h(X = 1) != 0: alg = 5.
 (4) The quotient factor (x^2 - |xi|^2/4)^2 does not vanish at r = |xi| = 1.
"""
import argparse, datetime, hashlib, json, os, platform, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
import sympy as sp  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); args = ap.parse_args(); log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    say('Exact algebraic multiplicity at the endpoint crossings')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__}')
    for f in ['audit/endpoint_multiplicity.py', 'data/m2_dispersion.json']:
        say(f'   {hashlib.sha256(open(os.path.join(ROOT, f), "rb").read()).hexdigest()}  {f}')
    x, s, e, b, X, v = sp.symbols('x s epsilon b X v')
    d = json.load(open(os.path.join(ROOT, 'data/m2_dispersion.json')))
    terms = [(int(i), int(j), int(k), sp.Rational(c)) for i, j, k, c in d['Q']]
    rec('(1) Q has total degree <= 4 in (x, s)', all(i + j <= 4 for i, j, k, c in terms))
    Q = sum(c * x ** i * s ** j * e ** k for i, j, k, c in terms)
    Qh = sp.expand(sum(c * x ** i * s ** j * b ** (4 - i - j) * e ** k for i, j, k, c in terms))
    rec('    Q_h(x, s, 1) == Q', sp.expand(Qh.subs(b, 1) - Q) == 0)
    rec('    lc_x(Q_h) == lc_x(Q) (x^4 coefficient free of s and b)',
        sp.expand(sp.Poly(Qh, x).LC() - sp.Poly(Q, x).LC()) == 0 and not sp.Poly(Q, x).LC().has(s))
    VL = (sp.Integer(34404) - 12 * sp.sqrt(5659689)) / 823; VU = (sp.Integer(58228) - 20 * sp.sqrt(1645345)) / 823
    for name, poly, vv in [('(2) lower endpoint, xi = (0,0,1,0), eps = eps_L', Q.subs(s, 0), VL),
                           ('(3) upper endpoint, xi = (1,0,0,0), eps = eps_U', Qh.subs({s: 1, b: 0}), VU)]:
        PX = sp.expand(poly.subs(x ** 4, X ** 2).subs(x ** 2, X))
        assert not PX.has(x)
        PXv = sp.expand(PX.subs(e ** 2, v)); assert not PXv.has(e)
        val = sp.radsimp(sp.expand(PXv.subs({X: 1, v: vv})))
        der = sp.radsimp(sp.expand(sp.diff(PXv, X).subs({X: 1, v: vv})))
        lc = sp.radsimp(sp.expand(sp.Poly(PXv, X).LC().subs(v, vv)))
        say(f'   {name}: Q(X=1) = {sp.simplify(val)}, dQ/dX(X=1) = {sp.nsimplify(der)} ~ {sp.N(der, 8)}, lc = {sp.N(lc, 8)}')
        rec(f'{name}: X = 1 is a simple root (value 0, derivative != 0, lc != 0)', sp.simplify(val) == 0 and sp.simplify(der) != 0 and sp.simplify(lc) != 0)
    rec('(4) quotient factor at r = 1: 1 - 1/4 != 0', sp.Rational(3, 4) != 0)
    say('   -> at both endpoint crossings, r = 1 has algebraic multiplicity exactly 5 on ker L.')
    line = 'ENDPOINT MULTIPLICITY: ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
