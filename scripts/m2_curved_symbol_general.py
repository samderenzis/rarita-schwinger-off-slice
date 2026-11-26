"""m2_curved_symbol_general.py -- curved-symbol identification for EVERY covector and EVERY charge (Milestone 2 input).

Usage:  python3 -I m2_curved_symbol_general.py [--out LOG]

As parent check_5 / step5_curved_symbol.py, but with the chart covector k = (k1, k2, k3, k4) fully symbolic and
q = q_slice + eta/35736 with eta symbolic (witness m, lam, b1, b2, pi' = 1/2). At y0 the principal symbol of the curved
eliminated operator K' equals i P(xi), xi = E(y0)^T k = (k1, k2, k3/3, k4/3), as an identity of rational functions
in (k, eta). Exact; SymPy.
"""
import argparse, datetime, hashlib, os, platform, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
from sympy import I, Matrix  # noqa: E402
import rs_model as rm  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out'); args = ap.parse_args(); log = []; T0 = time.time()

    def say(s=''):
        print(s, flush=True); log.append(s)
    say('Milestone 2 input -- curved-symbol identification with symbolic covector k and symbolic charge')
    say(f'date (UTC): {datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")}')
    say(f'python {platform.python_version()} | sympy {sp.__version__}')
    say('command: python3 -I scripts/m2_curved_symbol_general.py')
    for f in ['rs_model.py', 'm2_curved_symbol_general.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    eta = sp.Symbol('eta'); kk = list(sp.symbols('k1:5', real=True)); Lam = sp.Symbol('Lam', real=True)
    W0 = rm.witness_model(); qsl = rm.pauli_q(W0.m, W0.lam); W = rm.witness_model(q=qsl + eta / 35736)
    xi = [kk[0] * W.b1, kk[1] * W.b1, kk[2] * W.b2, kk[3] * W.b2]      # sin th_i = 1 at y0
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
    say(f'   {"OK  " if deg_ok else "FAIL"} degree in Lam <= 1: {deg_ok}')
    say(f'   {"OK  " if ident else "FAIL"} principal symbol of curved K\' at y0 == i P(E(y0)^T k) for symbolic k and eta: {ident}')
    say(f'   ({time.time() - T0:.1f} s)')
    ok = deg_ok and ident
    say('GENERAL CURVED-SYMBOL IDENTIFICATION: ' + ('PASS' if ok else 'FAIL'))
    if args.out:
        open(args.out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
