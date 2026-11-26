"""points.py -- exact definitions of the off-slice test points (pure Python, fractions only; no sympy).

Common parameters: the parent witness values lam^2 = 1/2, m^2 = 823/142944, b1 = 1, b2 = 1/3, pi' = 1/2, except the charge
    q = (1 + eps) * q_slice,   q_slice = lam*sqrt(2/3)*m = sqrt(2450894)/35736,   eps = SIGN * sqrt(V).
Equivalently q - q_slice = ETA/35736 with ETA = eps*sqrt(2450894) = ETA_RAT * sqrt(DELTA) (ETA_RAT carries the sign),
DELTA square-free. The symbol P depends on q only through ETA (step 4, part A), so every matrix at a point has entries in
    K = Q(i, sqrt(DELTA)),  basis 1, i, s, i*s  with s^2 = DELTA.
Covector (frame components at y0): XI = (a1, a2, b1, b2), all rational, with a1^2 + a2^2 = U (1 + U), b1^2 + b2^2 = 1 + U,
so U = |xi_a|^2/|xi_b|^2 and the light-branch eigenvalue R = |XI| = 1 + U is rational.

How the points were chosen (step 4, notes/STEP4_EXPLORATION.md): U is a rational point near u* = 400/9801 on the conic
w^2 = 4467 (u+1)(82875u+1267), with U and 1 + U sums of two rational squares; then V = eps^2 is the small positive root of
the crossing polynomial C(U, eps^2) = 0, and XI comes from two-square decompositions. For given eps^2 the crossing
polynomial has two roots u_- < u_+ (the two crossings into which the parent crossing splits); BRANCH says which one U is.
  A  : lower crossing, eps > 0          A- : same U, XI, V as A, eps < 0 (charge below the slice)
  B  : upper crossing, eps > 0          B- : same U, XI, V as B, eps < 0
None of this derivation is needed to CHECK a point: the step 5 scripts verify the crossing directly.
"""
from fractions import Fraction as Fr

D_SLICE = 2450894                         # q_slice = sqrt(D_SLICE)/35736
B1, B2 = Fr(1), Fr(1, 3)                  # sphere scales: S^2(1) x S^2(3)
PIP = Fr(1, 2)
LAM2, M2 = Fr(1, 2), Fr(823, 142944)

_A = dict(DELTA=50234001393, V=Fr(81608, 27765334551), ETA_RAT=Fr(404, 33736737),
          U=Fr(13484820232, 330742630449), BRANCH='lower',
          XI=(Fr(42108589454171386, 227779105414740649), Fr(186294557517875126, 2050011948732665841),
              Fr(501338840, 1431786279), Fr(1371952277, 1431786279)),
          R=Fr(344227450681, 330742630449))
_B = dict(DELTA=1589220303169, V=Fr(81608, 878393760583), ETA_RAT=Fr(404, 1067307121),
          U=Fr(426978392845, 10460177192916), BRANCH='upper',
          XI=(Fr(62858493798023514479, 310102446766203135396), Fr(52164981143920718591, 1395461010447914109282),
              Fr(53883658045, 52829177742), Fr(591024931, 26414588871)),
          R=Fr(10887155585761, 10460177192916))
POINTS = {
    'A': dict(_A, SIGN=+1),
    'A-': dict(_A, SIGN=-1, ETA_RAT=-_A['ETA_RAT']),
    'B': dict(_B, SIGN=+1),
    'B-': dict(_B, SIGN=-1, ETA_RAT=-_B['ETA_RAT']),
}


def crossing_poly_coeffs(v):
    """Coefficients (u^2, u, 1) of the crossing polynomial C(u, eps^2 = v) from step 4 (exact)."""
    return (718347204 * (677329 * v ** 2 - 95843288 * v + 2732361984),
            159632712 * (677329 * v ** 2 - 76236136 * v - 1003622400),
            8868484 * (677329 * v ** 2 - 56628984 * v + 368640000))


def self_check(name):
    """Exact consistency of the constants of one point (fractions only)."""
    p = POINTS[name]; a1, a2, b1, b2 = p['XI']; U, V = p['U'], p['V']
    c2, c1, c0 = crossing_poly_coeffs(V)
    # branch: U is the smaller (lower) or larger (upper) root of c2 u^2 + c1 u + c0; the other root is c0/(c2 U)
    other = Fr(c0) / (c2 * U)
    return {
        'ETA^2 == D_SLICE * V  (ETA = eps*sqrt(D_SLICE))': p['ETA_RAT'] ** 2 * p['DELTA'] == D_SLICE * V,
        'sign of ETA == SIGN': (p['ETA_RAT'] > 0) == (p['SIGN'] > 0),
        'a1^2 + a2^2 == U*(1+U)': a1 ** 2 + a2 ** 2 == U * (1 + U),
        'b1^2 + b2^2 == 1+U': b1 ** 2 + b2 ** 2 == 1 + U,
        '|XI|^2 == R^2, R == 1+U': sum(x * x for x in p['XI']) == p['R'] ** 2 and p['R'] == 1 + U,
        'V > 0 (eps real)': V > 0,
        'C(U, V) == 0 (U is a crossing for eps^2 = V; from the step 4 formula)': c2 * U * U + c1 * U + c0 == 0,
        f'U is the {p["BRANCH"]} root': (U < other) == (p['BRANCH'] == 'lower'),
        'q_slice^2 == lam^2 (2/3) m^2  (q_slice = sqrt(2450894)/35736)': Fr(D_SLICE, 35736 ** 2) == LAM2 * Fr(2, 3) * M2,
    }


if __name__ == '__main__':
    import sys
    for name in (sys.argv[1:] or POINTS):
        print(f'point {name}:')
        for k, v in self_check(name).items():
            print(f'   {k}: {v}')
