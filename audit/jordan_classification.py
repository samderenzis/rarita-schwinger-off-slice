"""jordan_classification.py -- exact Jordan classification at the certified off-slice crossings (SymPy, exact).

Usage:  python3 -I audit/jordan_classification.py [--out LOG]

Along xi = (s, 0, 1, 0), u = s^2, v = eps^2 (witness m, lam, b; q = (1 + eps) q_slice), with Q from
data/m2_dispersion.json (written as a quadratic in X = x^2: Q = a X^2 + b X + c):
  * the crossing condition is C(u, v) := Q(X = 1 + u) = 0;
  * at a crossing, r = sqrt(1 + u) is a root of Q(x) of multiplicity 1 iff disc_X Q (u, v) != 0 (then X = r^2 is a
    simple root of the quadratic and r != 0); multiplicity 2 iff disc_X Q = 0 (both sheets meet the light cone).
This script proves disc_X Q != 0 at every crossing of the certified ranges:
  lower branch 0 < v <= (266/100)^2, upper branch 0 < v <= (628/100)^2.
  (1) Any common zero of disc and C is a root of R(v) = Res_u(disc, C). R = const * v^2 * p5(v)^2 with p5 irreducible
      of degree 5. v > 0 excludes the factor v^2.
  (2) Exact real-root isolation of p5: its smallest positive root v0 (~ 32.2917, eps ~ 5.6826) lies below
      (628/100)^2 = 39.4384 and above (266/100)^2; all others lie above (628/100)^2. So on the lower range disc != 0 at
      every crossing; on the upper range only v in the isolating interval I0 of v0 needs care.
  (3) For v in I0 = [va, vb]: the upper root u_+(v) lies in [ulo, uhi] (exact Sturm count at va, no roots on the sides
      for v in I0, branch = larger root), and disc has no zero on [ulo, uhi] x [va, vb] (exact rational interval bound).
Conclusion: at every certified crossing the algebraic multiplicity of r on ker L is EXACTLY 5 (light factor 4,
Q factor 1; the quotient factor (x^2 - |xi|^2/4)^2 does not vanish at r). With rank (P - r)|ker L = 8 (>= 8 by the
box certificates, <= 8 identically by audit/restricted_symbol_fractions.py) the geometric multiplicity is exactly 4,
so the Jordan form of P|ker L on the generalised r-eigenspace is J2 + J1 + J1 + J1.
Also reported (exploratory context, outside the certified ranges): where disc and C have common real zeros.
"""
import argparse, datetime, hashlib, json, os, platform, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, '..')
import sympy as sp  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); args = ap.parse_args(); log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    say('Exact Jordan classification at the certified crossings')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__}')
    for f in ['audit/jordan_classification.py', 'data/m2_dispersion.json']:
        say(f'   {hashlib.sha256(open(os.path.join(ROOT, f), "rb").read()).hexdigest()}  {f}')
    x, s, e, X, u, v = sp.symbols('x s epsilon X u v')
    d = json.load(open(os.path.join(ROOT, 'data/m2_dispersion.json')))
    Q = sum(sp.Rational(c) * x ** int(i) * s ** int(j) * e ** int(k) for i, j, k, c in d['Q'])
    QX = sp.expand(Q.subs(x ** 4, X ** 2).subs(x ** 2, X)); rec('Q is a polynomial in X = x^2', not QX.has(x))
    QXv = sp.expand(QX.subs(s ** 2, u).subs(e ** 2, v)); rec('Q is a polynomial in u = s^2 and v = eps^2', not QXv.has(s) and not QXv.has(e))
    a, b, c = [sp.Poly(QXv, X).coeff_monomial(X ** k) for k in (2, 1, 0)]
    disc = sp.expand(b ** 2 - 4 * a * c); C = sp.expand(QXv.subs(X, 1 + u))
    say(f'   disc_X Q = {sp.factor(disc)}')
    say(f'   C(u, v) = {sp.factor(C)}')
    R = sp.factor(sp.resultant(disc, C, u))
    p5 = sp.Poly(377571474600343 * v ** 5 - 141153925904173916 * v ** 4 + 20332537164974086802 * v ** 3
                 - 1400990293471261118908 * v ** 2 + 45801950198017653991053 * v - 562553552914350795281250, v)
    const = sp.Poly(R, v).LC() / p5.LC() ** 2
    rec('Res_u(disc, C) == const * v^2 * p5(v)^2', sp.expand(R - const * v ** 2 * p5.as_expr() ** 2) == 0)
    rec('p5 irreducible over Q (SymPy)', p5.is_irreducible)
    c5 = [int(z) for z in p5.all_coeffs()]
    f11 = sp.Poly(p5.as_expr(), v, modulus=11).factor_list()[1]
    rec('    certificate: content 1, LC = 823^5 not divisible by 11, and p5 mod 11 is irreducible of degree 5',
        sp.gcd_list(c5) == 1 and c5[0] % 11 != 0 and len(f11) == 1 and f11[0][1] == 1 and f11[0][0].degree() == 5)
    lo2, up2 = Fr(266, 100) ** 2, Fr(628, 100) ** 2
    rts = p5.intervals()          # exact isolating intervals
    say(f'   isolating intervals of the real roots of p5: {[(str(I[0][0]), str(I[0][1])) for I in rts]}')
    below_lo = [I for I in rts if I[0][0] <= sp.Rational(lo2.numerator, lo2.denominator)]
    in_up = [I for I in rts if I[0][0] < sp.Rational(up2.numerator, up2.denominator)]
    rec('no root of p5 in (0, (2.66)^2]', len(below_lo) == 0 and all(I[0][0] > 0 for I in rts))
    rec('exactly one root v0 of p5 below (6.28)^2, with isolating interval inside ((2.66)^2, (6.28)^2)',
        len(in_up) == 1 and in_up[0][0][1] < sp.Rational(up2.numerator, up2.denominator) and in_up[0][0][0] > sp.Rational(lo2.numerator, lo2.denominator))
    va, vb = in_up[0][0]
    # refine the interval
    I0 = p5.refine_root(va, vb, eps=sp.Rational(1, 10 ** 12)); va, vb = I0
    say(f'   v0 in [{va}, {vb}]  (eps0 ~ {sp.N(sp.sqrt(va), 12)})')
    Cu = lambda vv: sp.Poly(C.subs(v, vv), u)
    ulo, uhi = sp.Rational(17, 10), sp.Rational(19, 10)
    n_in = Cu(va).count_roots(ulo, uhi); n_above = Cu(va).count_roots(uhi, 10 ** 9)
    side = all(sp.Poly(C.subs(u, ub), v).count_roots(va, vb) == 0 for ub in (ulo, uhi))
    rec('upper crossing u_+(v) in [17/10, 19/10] for all v in I0 (one root at va, none above, none on the sides)',
        n_in == 1 and n_above == 0 and side)
    n_neg = Cu(va).count_roots(sp.Rational(-2, 100), sp.Rational(-1, 100)); n_tot = Cu(va).count_roots()
    side2 = all(sp.Poly(C.subs(u, ub), v).count_roots(va, vb) == 0 for ub in (sp.Rational(-2, 100), sp.Rational(-1, 100)))
    rec('on I0 the only other real root of C(., v) is the continued lower branch, in [-2/100, -1/100] (u < 0: no real covector)',
        n_neg == 1 and n_tot == 2 and side2)
    # interval bound of disc on the box (u, v > 0): monomials monotone
    P = sp.Poly(disc, u, v); lo = hi = sp.Integer(0)
    for (i, j), cf in P.terms():
        mlo, mhi = ulo ** i * va ** j, uhi ** i * vb ** j
        lo += cf * (mlo if cf > 0 else mhi); hi += cf * (mhi if cf > 0 else mlo)
    rec(f'disc_X Q has no zero on [17/10, 19/10] x I0 (enclosure [{sp.N(lo, 6)}, {sp.N(hi, 6)}] excludes 0)', lo > 0 or hi < 0)
    VL = (sp.Integer(34404) - 12 * sp.sqrt(5659689)) / 823; VU = (sp.Integer(58228) - 20 * sp.sqrt(1645345)) / 823
    rec('no root of p5 in (0, eps_L^2] and the only root in (0, eps_U^2] is v0 (eps_L^2 < 32 <= v0 interval; eps_U^2 < 54 <= next root)',
        sp.N(VL) < 32 and (VL - 32).is_negative and (VU - 54).is_negative and rts[1][0][0] >= 54)
    say('   -> disc_X Q != 0 at every crossing with 0 < v <= eps_L^2 (lower) or 0 < v <= eps_U^2 (upper), endpoints included.')
    say('   -> algebraic multiplicity of r on ker L is exactly 5 there; with geometric multiplicity 4: J2 + 3 J1.')
    say('Exploratory context (outside the certified ranges): common real zeros of disc and C')
    for k in range(5):
        r_ = sp.CRootOf(p5.as_expr(), k)
        if not r_.is_real:
            continue
        vv = sp.N(r_, 40)
        cu = sp.Poly(C.subs(v, vv), u).nroots(n=30); du = sp.Poly(disc.subs(v, vv), u).nroots(n=30)
        com = [z for z in cu if min(abs(z - w) for w in du) < 1e-15]
        say(f'   eps ~ {sp.N(sp.sqrt(vv), 10)}: crossing roots u {[sp.N(z, 8) for z in cu]}; common with disc {[sp.N(z, 8) for z in com]}'
            f'{"  (real covector: triple crossing)" if com and sp.re(com[0]) > 0 else "  (u < 0: no real covector)"}')
    line = 'JORDAN CLASSIFICATION (exact): ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
