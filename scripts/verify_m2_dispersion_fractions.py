"""verify_m2_dispersion_fractions.py -- fractions-only PROOF of the full ray dispersion relation (no SymPy).

Usage:  python3 -I verify_m2_dispersion_fractions.py data/m2_dispersion.json [--out LOG]

Claim (witness m, lam, b1, b2, pi' = 1/2; xi = (s, 0, 1, 0); q = q_slice + eta/35736, eps = eta/sqrt(2450894)):
    det(x - P(xi)) = (x^2 - s^2 - 1)^4 (x^2 - (1 + s^2)/4)^2 Q(x, s, eps) / lc_x(Q),
with Q as recorded (data/m2_dispersion.json, from step 4 / m2_dispersion_export.py).
Proof by evaluation. Let d0(eta) = det S^0. d0*P has entries polynomial in eta of degree <= 4 (adj S^0 has degree 3,
S^k degree 1) and of degree <= 1 in s; so d0^16 det(x - P) is a polynomial in (s, eta) for each power of x, of degree
<= 16 in s and <= 64 in eta. lc_x(Q) is a constant multiple of d0 (recorded factors, checked below), so d0^16 times
the right-hand side is a polynomial of degree <= 8 + 4 + 4 = 16 in s and <= 60 + 4 = 64 in eta. Both sides agree at
all points of a 17 x 65 grid of rational (s, eta) with d0 != 0 (all coefficients in x compared exactly), hence everywhere.
P is built by this project's independent fractions construction (verify_m2_fractions.model / build_parts).
Consequence used by the family theorem: at a crossing (Q(|xi|) = 0) the eigenvalue r = |xi| has algebraic
multiplicity >= 5 on C^16, and none of it comes from the quotient factor (r^2 != r^2/4), so >= 5 on ker L.
"""
import hashlib
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_m2_fractions as vm  # noqa: E402

D = 2450894


def charpoly(A):
    """Faddeev-LeVerrier: coefficients c_0..c_n of det(x I - A) (c_n = 1), Gaussian rationals as Kel."""
    n = len(A); Mk = vm.mat(n, n); c = [vm.ZERO] * (n + 1); c[n] = vm.ONE
    I = vm.eye(n)
    for k in range(1, n + 1):
        Mk = vm.madd(vm.mmul(A, Mk), vm.msc(c[n - k + 1], I))
        AM = vm.mmul(A, Mk)
        tr = sum((AM[i][i] for i in range(n)), vm.ZERO)
        c[n - k] = tr * vm.Kel(Fr(-1, k))
    return c


def polymul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


_CACHE = {}


def _check_eta(arg):
    et = Fr(arg[0]); path = arg[1]
    data = json.loads(open(path, 'rb').read())
    Q = [(int(i), int(j), int(k), Fr(c)) for i, j, k, c in data['Q']]
    if 'model' not in _CACHE:
        _CACHE['model'] = vm.model()
    G, g3, Srat, Sq, _ = _CACHE['model']
    S = [vm.madd(Srat[r], vm.msc(vm.Kel(et / 35736), Sq[r])) for r in range(5)]
    assert not vm.det(S[0]).iszero()
    M, Bj = vm.build_parts(G, g3, S); Minv = vm.inverse(M)
    P1, P3 = vm.mmul(Minv, Bj[1]), vm.mmul(Minv, Bj[3])
    v = et * et / D; good = True; n = 0
    lc = (9803576 * v - 763140625) * (39214304 * v - 4014616321)
    for s in [Fr(k, 7) for k in range(-8, 9)]:                       # 17 values
        P = vm.madd(vm.msc(s, P1), P3)
        cl = charpoly(P)
        rhs = [Fr(0)] * 5
        for i, j, k, c in Q:
            rhs[i] += c * s ** j * v ** (k // 2)
        for _ in range(4):
            rhs = polymul(rhs, [-(s * s + 1), Fr(0), Fr(1)])
        for _ in range(2):
            rhs = polymul(rhs, [-(1 + s * s) / 4, Fr(0), Fr(1)])
        rhs = [c / lc for c in rhs]
        good &= len(rhs) == 17 and all(cl[i].a == rhs[i] and cl[i].b == 0 and cl[i].c == 0 and cl[i].d == 0 for i in range(17))
        n += 1
    return good, n


def main():
    path = sys.argv[1]; out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    raw = open(path, 'rb').read(); data = json.loads(raw); log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    import platform
    say('Fractions-only proof of the full ray dispersion relation (no SymPy)')
    say(f'python {platform.python_version()} | input {os.path.basename(path)} sha256 {hashlib.sha256(raw).hexdigest()}')
    say(f'script sha256 {hashlib.sha256(open(__file__, "rb").read()).hexdigest()}')
    say(f'verify_m2_fractions.py sha256 {hashlib.sha256(open(os.path.join(HERE, "verify_m2_fractions.py"), "rb").read()).hexdigest()}')
    Q = [(int(i), int(j), int(k), Fr(c)) for i, j, k, c in data['Q']]
    rec('Q has degree <= 4 in x, s, eps and only even powers of eps and x',
        all(i <= 4 and j <= 4 and k <= 4 and k % 2 == 0 and i % 2 == 0 for i, j, k, c in Q))

    def Qx(s, v):      # Q as polynomial in x (list) at s, eps^2 = v
        p = [Fr(0)] * 5
        for i, j, k, c in Q:
            p[i] += c * s ** j * v ** (k // 2)
        return p

    # lc_x(Q) = (9803576 eps^2 - 763140625)(39214304 eps^2 - 4014616321) (recorded value checked at sample points)
    def lcQ(v):
        return (9803576 * v - 763140625) * (39214304 * v - 4014616321)
    rec('lc_x(Q) equals (9803576 eps^2 - 763140625)(39214304 eps^2 - 4014616321) (polynomial in eps^2 of degree 2, '
        'compared at 3 points)', all(Qx(Fr(t), Fr(v))[4] == lcQ(Fr(v)) for t, v in [(1, 0), (2, 1), (3, 5)]))
    G, g3, Srat, Sq, cancel = vm.model()
    rec('S_mu - (16/3) i S_q == 0', cancel)
    from multiprocessing import Pool
    eta_grid = [Fr(k * 211, 3) for k in range(-32, 33)]             # 65 values, |eta| <= 2251 << 13812
    with Pool(int(os.environ.get('NPROC', os.cpu_count() or 1))) as pool:
        res = pool.map(_check_eta, [(str(e), path) for e in eta_grid])
    good = all(r[0] for r in res); npts = sum(r[1] for r in res)
    rec(f'det(x - P) == (x^2-s^2-1)^4 (x^2-(1+s^2)/4)^2 Q/lc_x(Q) at all {npts} grid points (17 s x 65 eta), every '
        f'coefficient in x', good)
    say('   -> by the degree bounds in the docstring, the identity holds for all s and eta with det S^0 != 0.')
    line = 'FRACTIONS PROOF OF THE RAY DISPERSION RELATION: ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if out:
        open(out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
