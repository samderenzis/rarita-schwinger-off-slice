"""m2_dispersion_export.py -- export the exact ray dispersion polynomial Q (step 4) for the fractions-only proof
verify_m2_dispersion_fractions.py of the FULL 16 x 16 characteristic polynomial along the ray.

Usage:  python3 -I m2_dispersion_export.py [--out LOG] [--json data/m2_dispersion.json]

Along xi = (s, 0, 1, 0), eps symbolic (q = (1 + eps) q_slice), witness m, lam, b, pi' = 1/2: Q is the factor of
charpoly(P|ker L) = (x^2 - s^2 - 1)^4 Q / lc_x(Q) (step 4 route: P restricted to ker L, fast). The full identity
    charpoly(P) = (x^2 - s^2 - 1)^4 (x^2 - (1 + s^2)/4)^2 Q(x, s, eps) / lc_x(Q)
(the quotient factor is paper Prop. 5.1(c)) is then PROVED by the fractions engine on a degree-sufficient grid, not here.
Writes Q's coefficients [i_x, j_s, k_eps, value].
"""
import argparse, datetime, hashlib, json, os, platform, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
import rs_model as rm  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); ap.add_argument('--json'); args = ap.parse_args()
    log = []; T0 = time.time()

    def say(s=''):
        print(s, flush=True); log.append(s)
    say('Ray dispersion relation, full 16 x 16 characteristic polynomial (SymPy)')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__}')
    say('command: python3 -I scripts/m2_dispersion_export.py')
    for f in ['rs_model.py', 'm2_dispersion_export.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    s, e = sp.symbols('s epsilon', real=True); x = sp.Symbol('x')
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam); W = rm.witness_model(q=(1 + e) * qsl)
    P = W.Psym([s, 0, 1, 0]).applyfunc(lambda z: sp.radsimp(sp.expand(z)))
    L = W.Lsym([s, 0, 1, 0]); Kb = sp.Matrix.hstack(*L.nullspace())
    A = ((Kb.H * Kb).applyfunc(sp.simplify).inv() * Kb.H * P * Kb).applyfunc(lambda z: sp.cancel(sp.radsimp(sp.expand(z))))
    cp = sp.cancel(A.charpoly(x).as_expr())
    n, d = sp.fraction(sp.together(cp)); n = sp.expand(sp.radsimp(sp.expand(n)))
    fl = sp.factor_list(n, x, s, e)
    Qs = [f for f, m in fl[1] if m == 1 and sp.degree(f, x) == 4]
    ok = len(Qs) == 1 and any(m == 4 for f, m in fl[1])
    Q = Qs[0]; lc = sp.Poly(Q, x).LC()
    ident = sp.cancel(cp - (x ** 2 - s ** 2 - 1) ** 4 * Q / lc) == 0
    say(f'   lc_x(Q) = {sp.factor(lc)}')
    say(f'   {"OK  " if ok and ident else "FAIL"} charpoly(P|ker L) == (x^2-s^2-1)^4 Q/lc_x(Q): {ok and ident}  [{time.time() - T0:.1f} s]')
    if args.json:
        pq = sp.Poly(Q, x, s, e)
        json.dump({'description': 'Q(x, s, eps): coefficients [i_x, j_s, k_eps, value]; eps = eta/sqrt(2450894).',
                   'Q': [[i, j, k, str(c)] for (i, j, k), c in pq.terms()]}, open(args.json, 'w'))
        say(f'   wrote {os.path.relpath(args.json)}  sha256 {hashlib.sha256(open(args.json, "rb").read()).hexdigest()}')
    say('DISPERSION EXPORT: ' + ('PASS' if ok and ident else 'FAIL'))
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok and ident else 1


if __name__ == '__main__':
    sys.exit(main())
