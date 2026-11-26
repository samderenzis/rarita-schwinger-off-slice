"""restricted_symbol_fractions.py -- exact checks behind RESTRICTED_SYMBOL_PROOF.md, in the fractions engine (no SymPy).

Usage:  python3 -I audit/restricted_symbol_fractions.py [--out LOG]

Witness m, lam, b1 = 1, b2 = 1/3, pi' = 1/2; q = q_slice + eta/35736 (eta rational here). P, L, S^rho built by the
independent fractions construction (scripts/verify_m2_fractions.py), never by SymPy or rs_model.

 (I)  Invariance of ker L:  L(xi) P(xi) = C1(xi) L(xi),  C1 = pi' sum_a xi_a Gamma^0 Gamma^a, for ALL xi in R^4 and
      all eta. Degree argument: d0(eta) = det S^0 clears every denominator of P (P = M^-1 B, B contains
      Gamma^{0ij} (S^0)^-1 S^k; adj S^0 has degree 3 and S^k degree 1 in eta), so d0 (L P - C1 L) has entries that are
      polynomials of degree <= 2 in each xi_j (L and P are linear in xi) and <= 4 in eta (d0 has degree 4, d0 P has
      degree <= 4). A grid with 3 values per xi_j and 5 values of eta (all with d0 != 0) therefore proves it.
 (II) C1(xi)^2 = pi'^2 |xi|^2 I and tr C1 = 0 for all xi (degree 2 in each xi_j: 3^4 grid), so the characteristic
      polynomial of C1 is (x^2 - pi'^2 |xi|^2)^2.
 (III) The ker L basis K(s) along the ray (data/m2_extended_certificate.json; entries a + b s) satisfies L K = 0 for
      all s and contains a 12 x 12 identity submatrix, so it has full column rank 12 = dim ker L (rank L = 4 by
      L L^dag = 3|xi|^2 I): it is a basis of ker L(s, 0, 1, 0) for every s.
 (IV) rank of M(tau, eta) = Dn (P - |xi|) K is <= 8 for ALL (tau, eta): every 9 x 9 minor is a polynomial of degree
      <= 9 * 4 = 36 in tau and in eta (each entry has degree <= 4 in each, proved in verify_m2_extended_fractions.py),
      so if rank <= 8 at all 37 x 37 points of a grid with distinct coordinates, every 9 x 9 minor vanishes on the grid
      and hence identically.
"""
import hashlib
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'scripts'))
import verify_m2_fractions as vm  # noqa: E402
import verify_point_fractions as vpf  # noqa: E402
import m2_box_tools as bt  # noqa: E402

ROOT = os.path.join(HERE, '..')


def main():
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    import platform
    say('Restricted-symbol checks (fractions engine, no SymPy)')
    say(f'python {platform.python_version()}')
    for f in ['audit/restricted_symbol_fractions.py', 'scripts/verify_m2_fractions.py', 'scripts/verify_point_fractions.py',
              'scripts/m2_box_tools.py', 'data/m2_extended_certificate.json']:
        say(f'   {hashlib.sha256(open(os.path.join(ROOT, f), "rb").read()).hexdigest()}  {f}')
    G, g3, Srat, Sq, cancel = vm.model()
    rec('S^rho = S_rat + (eta/35736) S_q exactly (the irrational mu cancels)', cancel)
    SP = [1, 2, 3, 4]; Z4 = vm.mat(4, 4); PIP = Fr(1, 2)

    def Lsym(xi):
        L = [[vm.ZERO] * 16 for _ in range(4)]
        for k in SP:
            blk = Z4
            for j in SP:
                blk = vm.madd(blk, vm.msc(xi[j - 1], g3(0, j, k)))
            for i in range(4):
                for c in range(4):
                    L[i][4 * (k - 1) + c] = blk[i][c]
        return L

    def C1(xi):
        return vm.msc(PIP, vm.madd(vm.madd(vm.msc(xi[0], vm.mmul(G[0], G[1])), vm.msc(xi[1], vm.mmul(G[0], G[2]))),
                              vm.madd(vm.msc(xi[2], vm.mmul(G[0], G[3])), vm.msc(xi[3], vm.mmul(G[0], G[4])))))

    say('(I) L P = C1 L for all xi, eta  (grid 3^4 x 5)')
    vals = [Fr(-1), Fr(0), Fr(2)]
    good = True; n = 0
    for et in [Fr(-4000), Fr(-1300, 7), Fr(0), Fr(2500, 3), Fr(9000)]:
        S = [vm.madd(Srat[r], vm.msc(vm.Kel(et / 35736), Sq[r])) for r in range(5)]
        assert not vm.det(S[0]).iszero()
        M, Bj = vm.build_parts(G, g3, S); Minv = vm.inverse(M)
        Pj = {j: vm.mmul(Minv, Bj[j]) for j in SP}
        for a in vals:
            for b in vals:
                for c in vals:
                    for d in vals:
                        xi = [a, b, c, d]
                        P = vm.mat(16, 16)
                        for j in SP:
                            if xi[j - 1] != 0:
                                P = vm.madd(P, vm.msc(xi[j - 1], Pj[j]))
                        L = Lsym(xi)
                        good &= vm.iszero_m(vm.madd(vm.mmul(L, P), vm.mmul(C1(xi), L), -vm.ONE)); n += 1
    rec(f'L P == C1 L at all {n} grid points (pi\' = 1/2)', good)
    say('(II) C1^2 = pi\'^2 |xi|^2 I, tr C1 = 0 (grid 3^4)')
    good = True
    for a in vals:
        for b in vals:
            for c in vals:
                for d in vals:
                    xi = [a, b, c, d]; C = C1(xi)
                    good &= vm.iszero_m(vm.madd(vm.mmul(C, C), vm.msc(PIP * PIP * sum(x * x for x in xi), vm.eye(4)), -vm.ONE))
                    good &= sum((C[i][i] for i in range(4)), vm.ZERO).iszero()
    rec('C1^2 == pi\'^2 |xi|^2 I and tr C1 == 0 at all 81 grid points', good)
    say('(III) K(s) has a 12 x 12 identity submatrix')
    data = json.load(open(os.path.join(ROOT, 'data/m2_extended_certificate.json')))
    Kc = data['K']
    rows_id = []
    for j in range(12):
        hit = [i for i in range(16) if Kc[i][j] == ['1', '0', '0', '0'] and all(
            Kc[i][jj] == ['0', '0', '0', '0'] for jj in range(12) if jj != j)]
        rows_id.append(hit[0] if hit else None)
    rec(f'rows {rows_id} of K form the 12 x 12 identity (constant in s)', None not in rows_id and len(set(rows_id)) == 12)
    good = True
    for sv in (Fr(0), Fr(1), Fr(-3, 2)):                         # L(s,0,1,0) K(s) has degree <= 2 in s
        Ks = [[vm.Kel(Fr(q_[0]), Fr(q_[1])) + vm.Kel(Fr(q_[2]), Fr(q_[3])) * sv for q_ in row] for row in Kc]
        good &= vm.iszero_m(vm.mmul(Lsym([sv, Fr(0), Fr(1), Fr(0)]), Ks))
    rec('L(s,0,1,0) K(s) == 0 for all s (degree <= 2 in s; checked at 3 values)', good)
    say('(IV) rank M(tau, eta) <= 8 identically (grid 37 x 37)')
    Mpoly = [[bt.poly_from_json(e) for e in row] for row in data['M']]
    taus = [Fr(k, 41) for k in range(1, 38)]                     # 37 distinct values in (0, 1)
    etas = [Fr(k * 400 - 7200, 3) for k in range(37)]            # 37 distinct values
    worst = 0; n = 0
    for t in taus:
        for e in etas:
            A = [[vm.Kel(*bt.poly_eval(Mpoly[i][j], t, e)) for j in range(12)] for i in range(16)]
            rk = vpf.rank(A); worst = max(worst, rk); n += 1
    rec(f'rank M <= 8 at all {n} grid points (max rank found {worst})', worst <= 8)
    say('   -> every 9 x 9 minor of M vanishes identically; with (III), rank (P - |xi|)|ker L <= 8 wherever Dn != 0.')
    line = 'RESTRICTED-SYMBOL CHECKS (fractions): ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if out:
        open(out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
