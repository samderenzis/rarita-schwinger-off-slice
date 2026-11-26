"""verify_point_fractions.py -- independent re-check of one off-slice point in pure Python (fractions only; no SymPy).

Usage:  python3 -I verify_point_fractions.py data/point_<NAME>_chain.json [--out LOG]

Second arithmetic engine for step 5. It does NOT use rs_model or SymPy. It rebuilds, from the paper's formulas:
  Gamma^A (the same 4x4 representation as the parent), F, X_A, S^rho, the paper's M = [Gamma^{0ik}], L(xi),
and works in the number field K = Q(i, s), s = sqrt(DELTA), elements a + b i + (c + d i) s with a, b, c, d in Q.

Irrational parameters. With m = sqrt(823/142944), l := lam*c1 = i*sqrt(3)/16 (lam = 1/sqrt2, c1 = i sqrt(3/128)),
  X_A = m Xm_A + l Xl_A,  Xm_A = -Gamma_A/3,  Xl_A = Fs Gamma_A - Gamma_A Fs/3,
  [X_a, X_b] = m^2 [Xm_a, Xm_b] + l^2 [Xl_a, Xl_b] + mu ([Xm_a, Xl_b] + [Xl_a, Xm_b]),  mu := m l,
with m^2 = 823/142944 and l^2 = -3/256 rational. The charge is q = q_slice + ETA/35736 with q_slice = -(16/3) i mu
(checked below: same square and same sign). So
  S^rho = S_rat + mu (S_mu - (16/3) i S_q) + (ETA/35736) S_q,   S_q := sum (i/2) Gamma^{ab rho} F_ab.
The script checks that S_mu - (16/3) i S_q = 0 exactly (the irrational mu cancels: this is where the slice
value of q enters), so S^rho has entries in K. Then, with B(xi) := [sum_j xi_j (Gamma^{ijk} - Gamma^{0ij} (S^0)^-1 S^k
+ pi' Gamma^0 Gamma^i Gamma^{0jk})]_{ik}, P = M^-1 B, and:
  checks: Clifford relations; det S^0 != 0 (and its value); L L^dag = 3|xi|^2 I; L v = L w = 0;
  (B - R M) v = 0 and (B - R M) w = M v  (i.e. (P - R) v = 0, (P - R) w = v); v != 0;
  dim ker (P - R)^j on ker L for j = 1, 2, 3, by Gaussian elimination over K, compared with the value recorded by
  step5_point_exact.py (the verifier does not assume an outcome). If no chain is recorded (NEGATIVE), only the ranks
  are checked.
"""
import hashlib
import json
import sys
from fractions import Fraction as Fr

DELTA = None  # set from the JSON


class Kel:
    """a + b i + (c + d i) s, s^2 = DELTA; a, b, c, d Fractions."""
    __slots__ = ('a', 'b', 'c', 'd')

    def __init__(self, a=0, b=0, c=0, d=0):
        self.a, self.b, self.c, self.d = Fr(a), Fr(b), Fr(c), Fr(d)

    def __add__(s, o):
        o = K(o); return Kel(s.a + o.a, s.b + o.b, s.c + o.c, s.d + o.d)
    __radd__ = __add__

    def __neg__(s):
        return Kel(-s.a, -s.b, -s.c, -s.d)

    def __sub__(s, o):
        return s + (-K(o))

    def __rsub__(s, o):
        return K(o) - s

    def __mul__(s, o):
        o = K(o)
        # (p1 + q1 s)(p2 + q2 s) = p1 p2 + DELTA q1 q2 + (p1 q2 + q1 p2) s, p, q Gaussian rationals
        def gm(x, y, u, v):
            return x * u - y * v, x * v + y * u
        p1p2 = gm(s.a, s.b, o.a, o.b); q1q2 = gm(s.c, s.d, o.c, o.d)
        p1q2 = gm(s.a, s.b, o.c, o.d); q1p2 = gm(s.c, s.d, o.a, o.b)
        return Kel(p1p2[0] + DELTA * q1q2[0], p1p2[1] + DELTA * q1q2[1], p1q2[0] + q1p2[0], p1q2[1] + q1p2[1])
    __rmul__ = __mul__

    def inv(s):
        # (p + q s)^-1 = (p - q s) / (p^2 - DELTA q^2), with the denominator in Q(i)
        n = s * Kel(s.a, s.b, -s.c, -s.d)
        assert n.c == 0 and n.d == 0
        den = n.a * n.a + n.b * n.b
        assert den != 0, 'division by zero in K'
        ninv = Kel(n.a / den, -n.b / den)
        return Kel(s.a, s.b, -s.c, -s.d) * ninv

    def __truediv__(s, o):
        return s * K(o).inv()

    def iszero(s):
        return s.a == 0 and s.b == 0 and s.c == 0 and s.d == 0

    def __eq__(s, o):
        return (s - K(o)).iszero()

    def conj(s):  # complex conjugation (s is real)
        return Kel(s.a, -s.b, s.c, -s.d)


def K(x):
    return x if isinstance(x, Kel) else Kel(x)


I_ = Kel(0, 1)
ZERO, ONE = Kel(), Kel(1)


# ---- matrices as lists of lists of Kel -------------------------------------------------------------------
def mat(n, m, f=lambda i, j: ZERO):
    return [[K(f(i, j)) for j in range(m)] for i in range(n)]


def eye(n):
    return mat(n, n, lambda i, j: ONE if i == j else ZERO)


def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum((A[i][t] * B[t][j] for t in range(k) if not A[i][t].iszero()), ZERO) for j in range(m)] for i in range(n)]


def madd(A, B, cb=ONE):
    return [[A[i][j] + cb * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def msc(c, A):
    return [[K(c) * x for x in row] for row in A]


def iszero_m(A):
    return all(x.iszero() for row in A for x in row)


def kron(A, B):
    return [[A[i // len(B)][j // len(B[0])] * B[i % len(B)][j % len(B[0])]
             for j in range(len(A[0]) * len(B[0]))] for i in range(len(A) * len(B))]


def H(A):
    return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]


def rref(A):
    """Reduced row echelon form over K; returns (R, pivots)."""
    R = [row[:] for row in A]; n, m = len(R), len(R[0]); piv = []; r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if not R[i][c].iszero()), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        iv = R[r][c].inv(); R[r] = [x * iv for x in R[r]]
        for i in range(n):
            if i != r and not R[i][c].iszero():
                f = R[i][c]; R[i] = [R[i][j] - f * R[r][j] for j in range(m)]
        piv.append(c); r += 1
        if r == n:
            break
    return R, piv


def rank(A):
    return len(rref(A)[1])


def inverse(A):
    n = len(A)
    R, piv = rref([A[i] + eye(n)[i] for i in range(n)])
    assert piv == list(range(n)), 'singular'
    return [row[n:] for row in R]


def nullspace(A):
    R, piv = rref(A); m = len(A[0]); free = [j for j in range(m) if j not in piv]; basis = []
    for f in free:
        v = [ZERO] * m; v[f] = ONE
        for r_, pc in enumerate(piv):
            v[pc] = -R[r_][f]
        basis.append(v)
    return [[basis[c][i] for c in range(len(basis))] for i in range(m)]   # columns


def det(A):
    A = [row[:] for row in A]; n = len(A); d = ONE
    for c in range(n):
        p = next((i for i in range(c, n) if not A[i][c].iszero()), None)
        if p is None:
            return ZERO
        if p != c:
            A[c], A[p] = A[p], A[c]; d = -d
        d = d * A[c][c]; iv = A[c][c].inv()
        for i in range(c + 1, n):
            if not A[i][c].iszero():
                f = A[i][c] * iv; A[i] = [A[i][j] - f * A[c][j] for j in range(n)]
    return d


def main():
    global DELTA
    path = sys.argv[1]; out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    raw = open(path, 'rb').read(); data = json.loads(raw)
    log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')

    import platform
    say(f'Independent fractions-only verification of off-slice point {data.get("point", "?")} (no SymPy)')
    say(f'python {platform.python_version()} | engine: fractions.Fraction + own arithmetic in Q(i, sqrt(DELTA))')
    say(f'input: {path.split("/")[-1]}  sha256 {hashlib.sha256(raw).hexdigest()}')
    say(f'script sha256 {hashlib.sha256(open(__file__, "rb").read()).hexdigest()}')
    DELTA = int(data['DELTA'])
    xi = [Fr(x) for x in data['XI']]; R = Fr(data['R']); den = Fr(data['DEN'])
    vec = lambda L: None if L is None else [[Kel(*[Fr(c) / den for c in e])] for e in L]
    v, w = vec(data['v']), vec(data['w'])
    ETA = Kel(0, 0, Fr(data['ETA_RAT']))   # ETA = ETA_RAT * s = eps*sqrt(2450894) = 35736 (q - q_slice)
    say(f'point {data["point"]}: {data["branch"]} crossing, sign of eps {data["sign_of_eps"]:+d}, eps^2 = {data["eps_squared"]}, DELTA = {DELTA}')
    PIP = Fr(1, 2); b1, b2 = Fr(1), Fr(1, 3)

    # Gamma matrices (same representation as the parent): G0 = i s3 x s0, G1 = s1 x s1, G2 = s1 x s2, G3 = s1 x s3, G4 = s2 x s0
    s0 = [[ONE, ZERO], [ZERO, ONE]]; s1 = [[ZERO, ONE], [ONE, ZERO]]
    s2 = [[ZERO, -I_], [I_, ZERO]]; s3 = [[ONE, ZERO], [ZERO, -ONE]]
    G = [msc(I_, kron(s3, s0)), kron(s1, s1), kron(s1, s2), kron(s1, s3), kron(s2, s0)]
    eta = [-1, 1, 1, 1, 1]; R5 = range(5); SP = [1, 2, 3, 4]; Z4 = mat(4, 4); E4 = eye(4)
    say('1. Clifford algebra')
    rec('{G^A, G^B} = 2 eta^{AB}', all(
        iszero_m(madd(madd(mmul(G[a], G[b]), mmul(G[b], G[a])), msc(-2 * (eta[a] if a == b else 0), E4)))
        for a in R5 for b in R5))

    def g3(a, b, c):
        return Z4 if len({a, b, c}) < 3 else mmul(mmul(G[a], G[b]), G[c])
    F = [[Fr(0)] * 5 for _ in R5]; F[1][2], F[2][1], F[3][4], F[4][3] = b1, -b1, b2, -b2
    Fs = Z4
    for a in R5:
        for b in R5:
            if F[a][b]:
                Fs = madd(Fs, msc(F[a][b], mmul(G[a], G[b])))
    Gl = [msc(eta[a], G[a]) for a in R5]
    Xm = [msc(Fr(-1, 3), Gl[a]) for a in R5]
    Xl = [madd(mmul(Fs, Gl[a]), msc(Fr(-1, 3), mmul(Gl[a], Fs))) for a in R5]
    m2, l2 = Fr(823, 142944), Fr(-3, 256)
    say('2. S^rho in K: the irrational mu = m*lam*c1 cancels against q_slice')
    # q_slice = -(16/3) i mu: compare squares (rational) and signs. mu = i*sqrt(823*3/(142944*256)), so
    # -(16/3) i mu = (16/3) sqrt(823*3/(142944*256)) > 0, and its square is (256/9)*823*3/(142944*256):
    rec('(16/3)^2 * 823*3/(142944*256) == q_slice^2 = lam^2 (2/3) m^2', Fr(256, 9) * Fr(823 * 3, 142944 * 256) == Fr(1, 2) * Fr(2, 3) * m2)
    Ein = [-(b1 ** 2 + b2 ** 2), -b2 ** 2, -b2 ** 2, -b1 ** 2, -b1 ** 2]

    def comm(A, B):
        return madd(mmul(A, B), mmul(B, A), -ONE)
    S, cancel_ok, Srat_l = [], True, []
    for r in R5:
        Srat = msc(Fr(1, 2) * Ein[r], G[r]); Smu = Z4; Sq = Z4
        for a in R5:
            for b in R5:
                g = g3(a, b, r)
                if iszero_m(g):
                    continue
                if F[a][b]:
                    Sq = madd(Sq, msc(Kel(0, Fr(1, 2)) * F[a][b], g))
                Srat = madd(Srat, msc(Fr(1, 2), mmul(g, madd(msc(m2, comm(Xm[a], Xm[b])), msc(l2, comm(Xl[a], Xl[b]))))))
                Smu = madd(Smu, msc(Fr(1, 2), mmul(g, madd(comm(Xm[a], Xl[b]), comm(Xl[a], Xm[b])))))
        cancel_ok &= iszero_m(madd(Smu, msc(Kel(0, Fr(-16, 3)), Sq)))
        S.append(madd(Srat, msc(ETA * Fr(1, 35736), Sq))); Srat_l.append(Srat)
    rec('S_mu - (16/3) i S_q == 0 for every rho (so S^rho has entries in K)', cancel_ok)
    def build(xi, S):
        S0i = inverse(S[0])
        M = [[ZERO] * 16 for _ in range(16)]; B = [[ZERO] * 16 for _ in range(16)]
        for i in SP:
            for k in SP:
                Mb = g3(0, i, k); Bb = Z4
                for j in SP:
                    term = madd(madd(g3(i, j, k), mmul(mmul(g3(0, i, j), S0i), S[k]), -ONE),
                                msc(PIP, mmul(mmul(G[0], G[i]), g3(0, j, k))))
                    Bb = madd(Bb, msc(xi[j - 1], term))
                for r_ in range(4):
                    for c in range(4):
                        M[4 * (i - 1) + r_][4 * (k - 1) + c] = Mb[r_][c]
                        B[4 * (i - 1) + r_][4 * (k - 1) + c] = Bb[r_][c]
        return M, B

    def Lsym(xi):
        L = [[ZERO] * 16 for _ in range(4)]
        for k in SP:
            blk = Z4
            for j in SP:
                blk = madd(blk, msc(xi[j - 1], g3(0, j, k)))
            for i in range(4):
                for c in range(4):
                    L[i][4 * (k - 1) + c] = blk[i][c]
        return L

    say('3. det S^0, constraint symbol')
    L = Lsym(xi)
    d0 = det(S[0])
    rec(f'det S^0 = {d0.a} (+ {d0.b} i + ({d0.c} + {d0.d} i) s), equal to the recorded value and != 0',
        not d0.iszero() and d0 == Kel(Fr(data['det_S0'])))
    rec('L L^dag == 3|xi|^2 I', iszero_m(madd(mmul(L, H(L)), msc(3 * sum(x * x for x in xi), E4), -ONE)))
    rec('|xi|^2 == R^2', sum(x * x for x in xi) == R * R)
    say('4a. Engine regression: the same code with ETA = 0 (Pauli slice) reproduces the parent witness')
    S_slice = [madd(Srat_l[r], Z4) for r in R5]
    d_sl = det(S_slice[0])
    rec('det S^0 on the slice == 3063716808343140625/132101812606622109696 (parent value)',
        d_sl == Kel(Fr(3063716808343140625, 132101812606622109696)))
    xs = [Fr(20, 99), Fr(0), Fr(1), Fr(0)]; rs = Fr(101, 99)
    pv = [11160500, 0, Kel(0, -16835940), Kel(0, -1495280), 0, 0, 0, 0, 55244475, 0, Kel(0, -83337903), Kel(0, -7401636), 0, 0, 0, 0]
    pw = [-Fr(3327596241179925, 5527168), -Fr(84691183525, 5568), Kel(0, Fr(407204312313361225, 679841664)),
          Kel(0, Fr(11902340395341825, 113306944)), Kel(0, -Fr(133046953528127925, 226613888)),
          Kel(0, -Fr(11666650486523175, 113306944)), -Fr(135032816027184075, 226613888), -Fr(1166781726825, 76096),
          Fr(552389505525, 19024), Fr(20361144375, 9512), 0, 0, 0, 0, 0, 0]
    pv = [[K(x)] for x in pv]; pw = [[K(x)] for x in pw]
    Ms, Bs = build(xs, S_slice); Ls = Lsym(xs); NBs = madd(Bs, msc(rs, Ms), -ONE)
    rec('parent chain: L v = L w = 0, (P-r) v = 0, (P-r) w = v (P built by this script)',
        iszero_m(mmul(Ls, pv)) and iszero_m(mmul(Ls, pw)) and iszero_m(mmul(NBs, pv))
        and iszero_m(madd(mmul(NBs, pw), mmul(Ms, pv), -ONE)))

    say('4b. Symbol and chain at the point')
    M, B = build(xi, S)
    rec('det M != 0', not det(M).iszero())
    NB = madd(B, msc(R, M), -ONE)          # M (P - R)
    if v is not None:
        rec('L v == 0', iszero_m(mmul(L, v)))
        rec('L w == 0', iszero_m(mmul(L, w)))
        rec('(P - R) v == 0   [M(P-R) v == 0]', iszero_m(mmul(NB, v)))
        rec('(P - R) w == v   [M(P-R) w == M v]', iszero_m(madd(mmul(NB, w), mmul(M, v), -ONE)))
        rec('v != 0', not iszero_m(v))
    else:
        say('   no chain recorded (NEGATIVE outcome): chain checks skipped')
    say('5. dim ker (P - R)^j restricted to ker L (Gaussian elimination over K)')
    kb = nullspace(L)
    rec('dim ker L == 12', len(kb[0]) == 12)
    Mi = inverse(M); N = mmul(Mi, NB)      # P - R
    dims = []; X = kb
    for j in range(1, 4):
        X = mmul(N, X); dims.append(12 - rank(X))
    say(f'   dims: {dims}')
    rec(f'dims == recorded {data["ker_dims_P_minus_R_on_kerL"][:3]}', dims == data['ker_dims_P_minus_R_on_kerL'][:3])
    pos = dims[1] > dims[0]
    say(f'   outcome from these ranks: {"POSITIVE (Jordan block in ker L)" if pos else "NEGATIVE (semisimple on ker L)"}')
    rec('outcome agrees with the recorded outcome', pos == data['outcome'].startswith('POSITIVE'))
    rec('chain recorded if and only if POSITIVE', (v is not None) == pos)
    line = f'FRACTIONS VERIFICATION OF POINT {data["point"]}: ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if out:
        open(out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
