"""witness_check.py -- the parent's witness_min.py checks, rebuilt on rs_model (q an independent parameter).

With q set back to the Pauli slice value it must print exactly the parent's golden witness_min.py
output (rarita-schwinger-well-posedness v1.0.11, results/golden_outputs.txt). This is a regression
test of the new code ON the slice; it says nothing yet about any point off the slice.

The chain v, w and the crossing (xi*, r) are copied verbatim from the parent's witness_min.py.

Usage:  python3 -I witness_check.py            (q passed explicitly as lam*sqrt(2/3)*m)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp  # noqa: E402
from sympy import I, Rational as Q, Matrix, eye  # noqa: E402
from rs_model import Model, G, E4, no_float  # noqa: E402

# parent witness parameters
lam2, m2, beta, pip = Q(1, 2), Q(823, 142944), Q(1, 3), Q(1, 2)
lam, m = sp.sqrt(lam2), sp.sqrt(m2)
q = lam * sp.sqrt(Q(2, 3)) * m                        # Pauli slice, passed EXPLICITLY
W = Model(m, lam, 1, beta, pip, q=q)
S, X = W.S, W.X

xi = [Q(20, 99), 0, 1, 0]; r = Q(101, 99)
P = W.Psym(xi); L = W.Lsym(xi); N = P - r * eye(16)

v = Matrix([11160500, 0, -16835940*I, -1495280*I, 0, 0, 0, 0, 55244475, 0, -83337903*I, -7401636*I, 0, 0, 0, 0])
w = Matrix([-Q(3327596241179925, 5527168), -Q(84691183525, 5568), Q(407204312313361225, 679841664)*I,
            Q(11902340395341825, 113306944)*I, -Q(133046953528127925, 226613888)*I, -Q(11666650486523175, 113306944)*I,
            -Q(135032816027184075, 226613888), -Q(1166781726825, 76096), Q(552389505525, 19024), Q(20361144375, 9512),
            0, 0, 0, 0, 0, 0])


def zero(Mx):
    return all(sp.expand(e) == 0 for e in Mx)


print("S^rho entries in Q(i) (sqrt(m), sqrt(lambda) cancelled):",
      all(sp.re(e).is_Rational and sp.im(e).is_Rational for Sx in S for e in Sx))
print("no Float atoms in G, X, S, P, L:", no_float(*(G + X + S + [P, L])))
print("det S0 =", sp.factor(S[0].det()))
print("L L^dag == 3|xi|^2 I:", zero(L * L.H - 3 * sum(x * x for x in xi) * E4))
print("L v == 0:", zero(L * v), "| L w == 0:", zero(L * w))
print("(P - r) v == 0:", zero(N * v), "| (P - r) w == v:", zero(N * w - v), "| v != 0:", not zero(v))
g = sp.symbols('gamma', positive=True); a = v - I * g * w
print("L a_gamma == 0:", zero(L * a), "| (gamma - i(P-r)) a_gamma == -i gamma^2 w:",
      zero((g * eye(16) - I * N) * a + I * g**2 * w))
print("sympy", sp.__version__)
