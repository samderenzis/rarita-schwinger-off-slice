"""verify_m2_fractions.py -- second engine for the Milestone 2 certificate (pure Python fractions; no SymPy).

Usage:  python3 -I verify_m2_fractions.py data/m2_minor.json [--out LOG]

Independent of SymPy and rs_model. Uses the arithmetic of verify_point_fractions.py (Gaussian rationals; DELTA unused
because eta is taken rational here) and rebuilds Gamma, S^rho(eta), M, B(xi) from the paper's formulas.
  1. Polynomial identity. For tau = 1/10 + a and eta = b on a 33 x 33 grid of rational points (a_k, b_l), it computes
     the 8 x 8 minor (recorded rows/cols) of Dn * (P(xi) - |xi|) K, xi = (2tau/(1-tau^2), 0, 1, 0), K as recorded, and
     compares with the recorded polynomial d(a, b). d has degree <= 32 in a and in b, and so does the computed minor
     (entries of Dn*(P-|xi|)K have degree <= 4 in tau and in eta after clearing, so an 8x8 minor has degree <= 32 in
     each: S^rho is affine in eta, det S^0 has degree 4 in eta and is cleared by Dn; P and K are affine in s and the
     (1-tau^2)^2 in Dn clears s^2 and s*|xi|), so agreement on the 33 x 33 tensor grid proves the identity.
  2. Taylor-box inequality: sum_{(i,j) != 0} (|Re c_ij| + |Im c_ij|) ha^i hb^j < |d(0,0)| for each recorded certificate.
  3. Box containment: on tau = 1/10 +- ha, the crossing polynomial (as a polynomial in v = eta^2/2450894) has no root
     v in [0, hb^2/2450894] (exact Sturm count).
"""
import hashlib
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_point_fractions as vpf  # noqa: E402

vpf.DELTA = 2          # irrelevant: all elements here have zero s-part
Kel, K_, ZERO, ONE, I_ = vpf.Kel, vpf.K, vpf.ZERO, vpf.ONE, vpf.I_
mmul, madd, msc, iszero_m, kron, inverse, det, mat, eye = (vpf.mmul, vpf.madd, vpf.msc, vpf.iszero_m, vpf.kron,
                                                            vpf.inverse, vpf.det, vpf.mat, vpf.eye)


def model():
    s0 = [[ONE, ZERO], [ZERO, ONE]]; s1 = [[ZERO, ONE], [ONE, ZERO]]
    s2 = [[ZERO, -I_], [I_, ZERO]]; s3 = [[ONE, ZERO], [ZERO, -ONE]]
    G = [msc(I_, kron(s3, s0)), kron(s1, s1), kron(s1, s2), kron(s1, s3), kron(s2, s0)]
    eta_s = [-1, 1, 1, 1, 1]; Z4 = mat(4, 4)

    def g3(a, b, c):
        return Z4 if len({a, b, c}) < 3 else mmul(mmul(G[a], G[b]), G[c])
    b1, b2 = Fr(1), Fr(1, 3)
    F = [[Fr(0)] * 5 for _ in range(5)]; F[1][2], F[2][1], F[3][4], F[4][3] = b1, -b1, b2, -b2
    Fs = Z4
    for a in range(5):
        for b in range(5):
            if F[a][b]:
                Fs = madd(Fs, msc(F[a][b], mmul(G[a], G[b])))
    Gl = [msc(eta_s[a], G[a]) for a in range(5)]
    Xm = [msc(Fr(-1, 3), Gl[a]) for a in range(5)]
    Xl = [madd(mmul(Fs, Gl[a]), msc(Fr(-1, 3), mmul(Gl[a], Fs))) for a in range(5)]
    m2, l2 = Fr(823, 142944), Fr(-3, 256)
    Ein = [-(b1 ** 2 + b2 ** 2), -b2 ** 2, -b2 ** 2, -b1 ** 2, -b1 ** 2]

    def comm(A, B):
        return madd(mmul(A, B), mmul(B, A), -ONE)
    Srat, Sq, cancel = [], [], True
    for r in range(5):
        sr = msc(Fr(1, 2) * Ein[r], G[r]); smu = Z4; sq = Z4
        for a in range(5):
            for b in range(5):
                g = g3(a, b, r)
                if iszero_m(g):
                    continue
                if F[a][b]:
                    sq = madd(sq, msc(Kel(0, Fr(1, 2)) * F[a][b], g))
                sr = madd(sr, msc(Fr(1, 2), mmul(g, madd(msc(m2, comm(Xm[a], Xm[b])), msc(l2, comm(Xl[a], Xl[b]))))))
                smu = madd(smu, msc(Fr(1, 2), mmul(g, madd(comm(Xm[a], Xl[b]), comm(Xl[a], Xm[b])))))
        cancel &= iszero_m(madd(smu, msc(Kel(0, Fr(-16, 3)), sq)))
        Srat.append(sr); Sq.append(sq)
    return G, g3, Srat, Sq, cancel


def build_parts(G, g3, S):
    """M (paper) and the xi-linear pieces of B: B(xi) = sum_j xi_j Bj."""
    Z4 = mat(4, 4); S0i = inverse(S[0]); SP = [1, 2, 3, 4]; PIP = Fr(1, 2)
    M = [[ZERO] * 16 for _ in range(16)]; Bj = {j: [[ZERO] * 16 for _ in range(16)] for j in SP}
    for i in SP:
        for k in SP:
            Mb = g3(0, i, k)
            for r_ in range(4):
                for c in range(4):
                    M[4 * (i - 1) + r_][4 * (k - 1) + c] = Mb[r_][c]
            for j in SP:
                term = madd(madd(g3(i, j, k), mmul(mmul(g3(0, i, j), S0i), S[k]), -ONE),
                            msc(PIP, mmul(mmul(G[0], G[i]), g3(0, j, k))))
                for r_ in range(4):
                    for c in range(4):
                        Bj[j][4 * (i - 1) + r_][4 * (k - 1) + c] = term[r_][c]
    return M, Bj


def sturm_count(p, lo, hi):
    """Number of distinct real roots of polynomial p (list of Fractions, low degree first) in (lo, hi]."""
    def trim(q):
        while q and q[-1] == 0:
            q = q[:-1]
        return q

    def ev(q, x):
        return sum(c * x ** i for i, c in enumerate(q))

    def rem(a, b):
        a = a[:]
        while len(a) >= len(b) and a:
            f = a[-1] / b[-1]; sh = len(a) - len(b)
            for i, c in enumerate(b):
                a[i + sh] -= f * c
            a = trim(a)
        return a
    p = trim(p); dp = trim([i * c for i, c in enumerate(p)][1:])
    seq = [p, dp]
    while seq[-1] and len(seq[-1]) > 1:
        r = rem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])

    def changes(x):
        vals = [ev(q, x) for q in seq if q]; vals = [v for v in vals if v != 0]
        return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))
    return changes(lo) - changes(hi)


def main():
    path = sys.argv[1]; out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    raw = open(path, 'rb').read(); data = json.loads(raw); log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    import platform
    say('Independent fractions-only verification of the Milestone 2 certificate (no SymPy)')
    say(f'python {platform.python_version()} | input {os.path.basename(path)} sha256 {hashlib.sha256(raw).hexdigest()}')
    say(f'script sha256 {hashlib.sha256(open(__file__, "rb").read()).hexdigest()}')
    rows, cols = data['rows'], data['cols']
    coeffs = {(int(i), int(j)): Kel(Fr(re), Fr(im)) for i, j, re, im in data['coefficients']}
    da = max(i for i, j in coeffs); db = max(j for i, j in coeffs)
    rec(f'recorded degrees (a, b) = ({da}, {db}) <= (32, 32)', da <= 32 and db <= 32)
    Kc = [[e for e in row] for row in data['K']]
    G, g3, Srat, Sq, cancel = model()
    rec('S_mu - (16/3) i S_q == 0 (S^rho = S_rat + (eta/35736) S_q exactly)', cancel)

    def d_poly(a, b):
        return sum((c * (a ** i) * (b ** j) for (i, j), c in coeffs.items()), ZERO)
    say('1. Polynomial identity on a 33 x 33 grid')
    grid_a = [Fr(k - 16, 400) for k in range(33)]          # a in [-0.04, 0.04]
    grid_b = [Fr(k - 16, 2) for k in range(33)]            # eta in [-8, 8]
    allgood = True
    for bb in grid_b:
        S = [madd(Srat[r], msc(Kel(bb / 35736), Sq[r])) for r in range(5)]
        M, Bj = build_parts(G, g3, S)
        Minv = inverse(M)
        for aa in grid_a:
            t = Fr(1, 10) + aa
            s = 2 * t / (1 - t * t); x = (1 + t * t) / (1 - t * t)
            Dn = (t * t - 1) ** 2 * (4 * bb * bb - 27625 ** 2) * (16 * bb * bb - 63361 ** 2)
            Bx = madd(msc(s, Bj[1]), Bj[3])                       # xi = (s, 0, 1, 0)
            # P - x = M^-1 B - x  ->  we need (P - x) K: compute M^-1 (B - x M) K
            Kt = [[Kel(Fr(e[0]), Fr(e[1])) + Kel(Fr(e[2]), Fr(e[3])) * s for e in row] for row in Kc]
            NB = madd(Bx, msc(x, M), -ONE)
            PK = mmul(Minv, mmul(NB, Kt))
            sub = [[PK[i][j] * Dn for j in cols] for i in rows]
            val = det(sub)
            if not (val == d_poly(aa, bb)):
                allgood = False
    rec('computed minor == recorded d(a, b) at all 1089 grid points (so identical as polynomials)', allgood)
    d00 = coeffs[(0, 0)]
    say('2. Taylor-box inequalities')
    for cert in data['certificates']:
        ha, hb = Fr(cert['ha']), Fr(cert['hb'])
        bnd = sum(((abs(c.a) + abs(c.b)) * ha ** i * hb ** j for (i, j), c in coeffs.items() if (i, j) != (0, 0)), Fr(0))
        rec(f'eps0 = {cert["eps0"]}: bound < |d(0,0)| (ratio {float(bnd / abs(d00.a)):.6f}); d(0,0) real', d00.b == 0 and bnd < abs(d00.a))
        # 3. containment: crossing polynomial in v on the box sides
        e0 = Fr(cert['eps0']); hb2 = hb * hb; D = 2450894
        rec(f'   hb^2 >= eps0^2 * 2450894', hb2 >= e0 * e0 * D)
        okside = True
        for t in (Fr(1, 10) - ha, Fr(1, 10) + ha):
            # C(u, v) = c2(v) u^2 + c1(v) u + c0(v), u = s^2 = 4t^2/(1-t^2)^2; times (1-t^2)^4:
            c2 = [Fr(718347204 * 2732361984), Fr(-718347204 * 95843288), Fr(718347204 * 677329)]
            c1 = [Fr(159632712 * -1003622400), Fr(-159632712 * 76236136), Fr(159632712 * 677329)]
            c0 = [Fr(8868484 * 368640000), Fr(-8868484 * 56628984), Fr(8868484 * 677329)]
            w = 1 - t * t
            g = [16 * t ** 4 * c2[k] + 4 * t * t * w * w * c1[k] + w ** 4 * c0[k] for k in range(3)]
            n = sturm_count(g, Fr(-1, 10 ** 30), hb2 / D)      # roots v in (-tiny, hb^2/D] covers [0, hb^2/D]
            okside &= n == 0
        rec(f'   crossing polynomial has no root on the box sides tau = 1/10 +- {ha} for |eta| <= hb', okside)
    line = 'FRACTIONS VERIFICATION OF MILESTONE 2 CERTIFICATE: ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if out:
        open(out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
