"""m2_rotation_invariance.py -- the symbols are equivariant under rotations within each sphere factor (Milestone 2 input).

Usage:  python3 -I m2_rotation_invariance.py [--out LOG]

For the rotation R by angle theta in the (1,2) plane (and separately the (3,4) plane) of the spatial frame, with
S = cos(theta/2) + sin(theta/2) Gamma^i Gamma^j and Shat = R^T (x) S acting on Phi = (psi_1..psi_4):
    S Gamma^a S^-1 = R^a_b Gamma^b (a = 1..4),  S Gamma^0 S^-1 = Gamma^0,
    P(R^T xi) Shat = Shat P(xi),   L(R^T xi) Shat = S L(xi),
as exact identities in symbolic xi = (xi1..xi4), symbolic charge (q = q_slice + eta/35736, eta symbolic), witness
m, lam, b1, b2, pi' = 1/2, and alpha = cos(theta/2), beta = sin(theta/2) symbolic modulo alpha^2 + beta^2 = 1.
Consequences: the spectrum of P(xi), ker L(xi), and the Jordan structure of P on ker L depend on xi only through
|xi_a|^2, |xi_b|^2 (and, by homogeneity P(c xi) = c P(xi), only through their ratio up to scaling). In particular every
xi with xi_b != 0 is similar to a multiple of a point (s, 0, 1, 0), s >= 0. Exact; SymPy.
"""
import argparse, datetime, hashlib, os, platform, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
import rs_model as rm  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); args = ap.parse_args(); log = []; ok_all = True; T0 = time.time()

    def say(s=''):
        print(s, flush=True); log.append(s)
    say('Milestone 2 input -- rotation equivariance of P and L')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__}')
    say('command: python3 -I scripts/m2_rotation_invariance.py')
    for f in ['rs_model.py', 'm2_rotation_invariance.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    eta, al, be = sp.symbols('eta alpha beta', real=True); xs = sp.symbols('xi1:5', real=True)
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam); W = rm.witness_model(q=qsl + eta / 35736)
    G = rm.G; E4 = sp.eye(4)

    def red(e):
        e = sp.expand(e)
        return sp.expand(sp.rem(sp.Poly(e, al), sp.Poly(al ** 2 + be ** 2 - 1, al)).as_expr()) if e.has(al) else e
    zero = lambda M: all(red(sp.numer(sp.together(sp.radsimp(sp.expand(z))))) == 0 for z in M)
    for (i, j) in [(1, 2), (3, 4)]:
        S = al * E4 + be * G[i] * G[j]; Si = al * E4 - be * G[i] * G[j]
        c, s = al ** 2 - be ** 2, 2 * al * be
        R = sp.eye(4); R[i - 1, i - 1] = c; R[i - 1, j - 1] = -s; R[j - 1, i - 1] = s; R[j - 1, j - 1] = c
        Sh = sp.kronecker_product(R.T, S); Shi = sp.kronecker_product(R, Si)
        xr = list(R.T * sp.Matrix(xs))
        checks = {
            'S S^-1 = 1': zero(S * Si - E4),
            'S Gamma^a S^-1 = R^a_b Gamma^b, S Gamma^0 S^-1 = Gamma^0':
                all(zero(S * G[a] * Si - sum((R[a - 1, bb - 1] * G[bb] for bb in rm.SP), sp.zeros(4))) for a in rm.SP)
                and zero(S * G[0] * Si - G[0]),
            'Shat Shat^-1 = 1': zero(Sh * Shi - sp.eye(16)),
            'R^T R = 1 (R a rotation)': zero(R.T * R - E4),
            'L(R^T xi) Shat == S L(xi)': zero(W.Lsym(xr) * Sh - S * W.Lsym(list(xs))),
            'P(R^T xi) Shat == Shat P(xi)': zero(W.Psym(xr) * Sh - Sh * W.Psym(list(xs))),
        }
        for k, v in checks.items():
            ok_all &= v; say(f'   {"OK  " if v else "FAIL"} rotation in the ({i},{j}) plane: {k}: {v}  [{time.time() - T0:.1f} s]')
    say('ROTATION EQUIVARIANCE: ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED'))
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
