"""m2_box_tools.py -- exact (fractions-only) tools for the extended Milestone 2 box certificate. No SymPy.

Objects
  * Gaussian rationals as pairs (re, im) of Fractions.
  * A bivariate polynomial p(tau, eta) as a dict {(i, j): (re, im)}.
  * M(tau, eta): the 16 x 12 polynomial matrix Dn * (P(xi) - |xi|) K, xi = (2tau/(1-tau^2), 0, 1, 0) (see
    m2_family_q_direction.py); its columns span (part of) ker L, so rank of (P - |xi|) on ker L >= rank M where Dn != 0.

Box certificate (one box B = [t_lo, t_hi] x [e_lo, e_hi] in (tau, eta)):
  rank-8 test: for an 8 x 8 submatrix A = M[rows, cols] and a rational matrix C, with the box centre (tc, ec) and
  half-widths (ha, hb): every entry of C A(tau, eta) - I is a polynomial; |entry| <= |value at centre| +
  sum_{(i,j) != 0} |shifted coeff_ij| ha^i hb^j (|z| <= |Re z| + |Im z|). If every row sum of these bounds is < 1,
  then ||C A - I||_inf < 1 on B, so C A and hence A are invertible on B: rank M >= 8 on B.
Crossing containment (one ε-segment, v = eta^2 = eps^2 * 2450894 in [v_lo, v_hi]):
  Ct(tau, v) = 16 tau^4 c2(v/D) + 4 tau^2 (1 - tau^2)^2 c1(v/D) + (1 - tau^2)^4 c0(v/D) is (1 - tau^2)^4 C(s^2, eps^2).
  If Ct(t_lo, v) and Ct(t_hi, v) have no root v in [v_lo, v_hi] and Ct(., v_lo) has exactly one root tau in
  (t_lo, t_hi), then the crossing branch through that root stays in (t_lo, t_hi) for all v in [v_lo, v_hi]
  (roots move continuously and cannot cross the sides). Branch identification: see segment_containment.
"""
from fractions import Fraction as Fr
from math import comb

D = 2450894


# ---------- Gaussian rationals --------------------------------------------------------------------------------
def gadd(x, y): return (x[0] + y[0], x[1] + y[1])
def gmul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def gabs_bound(x): return abs(x[0]) + abs(x[1])


# ---------- bivariate polynomials -----------------------------------------------------------------------------
def poly_from_json(terms):
    return {(int(i), int(j)): (Fr(re), Fr(im)) for i, j, re, im in terms}


def poly_eval(p, t, e):
    acc = (Fr(0), Fr(0))
    for (i, j), c in p.items():
        acc = gadd(acc, (c[0] * t ** i * e ** j, c[1] * t ** i * e ** j))
    return acc


def poly_shift(p, tc, ec):
    """Coefficients of p(tc + a, ec + b) in a, b."""
    out = {}
    for (i, j), c in p.items():
        for k in range(i + 1):
            fk = comb(i, k) * tc ** (i - k)
            for l in range(j + 1):
                f = fk * comb(j, l) * ec ** (j - l)
                o = out.get((k, l), (Fr(0), Fr(0)))
                out[(k, l)] = (o[0] + c[0] * f, o[1] + c[1] * f)
    return out


def box_rank8(Mpoly, rows, cols, C, tc, ec, ha, hb):
    """Return (ok, max row-sum bound of |C A - I|) on the box centred (tc, ec) with half-widths (ha, hb)."""
    sh = {(r, c): poly_shift(Mpoly[r][c], tc, ec) for r in rows for c in cols}
    worst = Fr(0)
    for i in range(8):
        rowsum = Fr(0)
        for j in range(8):
            ent = {}
            for k in range(8):
                for m, coef in sh[(rows[k], cols[j])].items():
                    o = ent.get(m, (Fr(0), Fr(0)))
                    ent[m] = gadd(o, gmul(C[i][k], coef))
            c00 = ent.get((0, 0), (Fr(0), Fr(0)))
            c00 = (c00[0] - (1 if i == j else 0), c00[1])
            bnd = gabs_bound(c00) + sum(gabs_bound(c) * ha ** a * hb ** b for (a, b), c in ent.items() if (a, b) != (0, 0))
            rowsum += bnd
        worst = max(worst, rowsum)
    return worst < 1, worst


# ---------- crossing polynomial and Sturm counts ---------------------------------------------------------------
def c_coeffs_in_v(e2):
    """(c2, c1, c0) of C(u, eps^2) at eps^2 = e2 (Fractions); from points.crossing_poly_coeffs (step 4)."""
    return (718347204 * (677329 * e2 ** 2 - 95843288 * e2 + 2732361984),
            159632712 * (677329 * e2 ** 2 - 76236136 * e2 - 1003622400),
            8868484 * (677329 * e2 ** 2 - 56628984 * e2 + 368640000))


def Ct_in_v(t, ray=1):
    """Ct(t, v) as a polynomial in v (list, low degree first), t a Fraction. ray = 2: second ray xi = (1, 0, t', 0),
    t' = 2t/(1 - t^2), where u' = t'^2 = 1/u and the roles of c2 and c0 are exchanged."""
    w = 1 - t * t
    a, b, c = 16 * t ** 4, 4 * t * t * w * w, w ** 4
    # c2 = 718347204*(677329 e^2 - 95843288 e + 2732361984), e = v/D, etc.
    def cv(k, p2, p1, p0):
        return [k * p0, k * Fr(p1, D), k * Fr(p2, D * D)]
    c2 = cv(718347204, 677329, -95843288, 2732361984)
    c1 = cv(159632712, 677329, -76236136, -1003622400)
    c0 = cv(8868484, 677329, -56628984, 368640000)
    if ray == 2:
        c2, c0 = c0, c2
    return [a * c2[n] + b * c1[n] + c * c0[n] for n in range(3)]


def Ct_in_tau(v, ray=1):
    """Ct(tau, v) as a polynomial in tau (list, low degree first, degree 8), v a Fraction (ray: see Ct_in_v)."""
    c2, c1, c0 = c_coeffs_in_v(v / D)
    if ray == 2:
        c2, c0 = c0, c2
    # 16 t^4 c2 + 4 t^2 (1 - t^2)^2 c1 + (1 - t^2)^4 c0
    p = [Fr(0)] * 9
    p[4] += 16 * c2
    for k, cf in enumerate([1, -2, 1]):          # (1 - t^2)^2 = 1 - 2t^2 + t^4
        p[2 + 2 * k] += 4 * cf * c1
    for k, cf in enumerate([1, -4, 6, -4, 1]):   # (1 - t^2)^4
        p[2 * k] += cf * c0
    return p


def _trim(q):
    q = list(q)
    while q and q[-1] == 0:
        q.pop()
    return q


def _ev(q, x):
    acc = Fr(0)
    for c in reversed(q):
        acc = acc * x + c
    return acc


def _rem(a, b):
    a = list(a)
    while len(a) >= len(b) and a:
        f = a[-1] / b[-1]; sh = len(a) - len(b)
        for i, c in enumerate(b):
            a[i + sh] -= f * c
        a = _trim(a)
    return a


def sturm_count(p, lo, hi):
    """Number of distinct real roots of p in (lo, hi] (p not identically zero)."""
    p = _trim(p); dp = _trim([i * c for i, c in enumerate(p)][1:])
    seq = [p, dp] if dp else [p]
    while len(seq) > 1 and len(seq[-1]) > 1:
        r = _rem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])

    def changes(x):
        vals = [v for v in (_ev(q, x) for q in seq) if v != 0]
        return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))
    return changes(lo) - changes(hi)


def no_root_closed(p, lo, hi):
    """True iff p has no real root in the closed interval [lo, hi]."""
    p = _trim(p)
    return _ev(p, lo) != 0 and sturm_count(p, lo, hi) == 0


def segment_containment(t_lo, t_hi, v_lo, v_hi, branch=None, ray=1, touch_zero=False):
    """Exact containment test for one eps-segment (see module docstring). With branch = 'lower' ('upper') it also
    checks that the root in the box is the smallest (largest) root of Ct(., v_lo) in (0, 1): the only roots there are
    tau_- < tau_+ (Ct is even and palindromic in tau, roots +-tau, +-1/tau), corresponding to u_- < u_+."""
    if touch_zero:
        # End box [0, t_hi]: the branch is the SMALLEST root of Ct(., v) in (0, 1) (lower branch on ray 1; upper
        # branch on ray 2). It starts in (0, t_hi) at v_lo and cannot reach t_hi; it may leave only through 0, which is
        # where the branch ends (u = 0 on ray 1 at eps_L; u' = 0 on ray 2 at eps_U). See EPSILON_ENDPOINT_CLOSURE.md.
        if not (t_lo == 0 and 0 < t_hi < 1):
            return False
        p = Ct_in_tau(v_lo, ray)
        start = sturm_count(p, Fr(0), t_hi) == 1 and _ev(p, t_hi) != 0 and _ev(p, Fr(0)) != 0
        sides = no_root_closed(Ct_in_v(t_hi, ray), v_lo, v_hi)
        return start and sides
    if not (0 < t_lo < t_hi < 1):
        return False
    p = Ct_in_tau(v_lo, ray)
    start = sturm_count(p, t_lo, t_hi) == 1 and _ev(p, t_hi) != 0
    sides = no_root_closed(Ct_in_v(t_lo, ray), v_lo, v_hi) and no_root_closed(Ct_in_v(t_hi, ray), v_lo, v_hi)
    ident = True
    smallest = (branch == 'lower') if ray == 1 else (branch == 'upper')
    largest = (branch == 'upper') if ray == 1 else (branch == 'lower')
    if smallest:
        ident = sturm_count(p, Fr(0), t_lo) == 0
    elif largest:
        ident = sturm_count(p, t_hi, Fr(1)) == 0 and _ev(p, Fr(1)) != 0
    return start and sides and ident
