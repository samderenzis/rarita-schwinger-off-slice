"""verify_m2_gaps_fractions.py -- second check of the endpoint-gap certificate (pure fractions; no SymPy).

Usage:  python3 -I verify_m2_gaps_fractions.py data/m2_gap_certificate.json [--out LOG]

  1. M1 recorded == the extended certificate's M (already proved there) -- here re-proved on a 5 x 5 grid from the
     independent fractions construction of P (degree <= 4 in tau and eta).
  2. M2 (ray 2, xi' = (1, 0, t, 0), t = 2sigma/(1 - sigma^2)) recomputed from the independent construction on a 7 x 5
     grid in (sigma, eta): entries have degree <= 6 in sigma and <= 4 in eta (checked on the recorded polynomials;
     true entries: P affine in t and K'(t)(1 + t^2) of degree <= 3 in t give, after Dn' = (sigma^2 - 1)(sigma^2 + 1)^2 ...,
     degree <= 6 in sigma; eta as in ray 1). K' is also checked: L(1,0,t,0) K'(t) = 0 and the 12 x 12 identity block.
  3. Every box: rank-8 Taylor test with the independent Taylor-shift implementation (verify_m2_extended_fractions);
     containment with m2_box_tools.segment_containment (shared code); eta enclosures; contiguity; end boxes touch 0;
     end eps exceed eps_L, eps_U exactly.
"""
import hashlib
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_m2_fractions as vm  # noqa: E402
import verify_m2_extended_fractions as vx  # noqa: E402
import m2_box_tools as bt  # noqa: E402

D = bt.D


def main():
    path = sys.argv[1]; out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    raw = open(path, 'rb').read(); data = json.loads(raw); log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    import platform
    say('Second check of the endpoint-gap certificate (fractions only, no SymPy)')
    say(f'python {platform.python_version()} | input {os.path.basename(path)} sha256 {hashlib.sha256(raw).hexdigest()}')
    for f in ['verify_m2_gaps_fractions.py', 'verify_m2_extended_fractions.py', 'verify_m2_fractions.py', 'm2_box_tools.py']:
        say(f'   {hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()}  scripts/{f}')
    M1 = [[bt.poly_from_json(e) for e in row] for row in data['M1']]
    M2 = [[bt.poly_from_json(e) for e in row] for row in data['M2']]
    G, g3, Srat, Sq, cancel = vm.model(); rec('S^rho = S_rat + (eta/35736) S_q exactly', cancel)
    Z4 = vm.mat(4, 4)

    def Lsym(xi):
        L = [[vm.ZERO] * 16 for _ in range(4)]
        for k in [1, 2, 3, 4]:
            blk = Z4
            for j in [1, 2, 3, 4]:
                blk = vm.madd(blk, vm.msc(xi[j - 1], g3(0, j, k)))
            for i in range(4):
                for c in range(4):
                    L[i][4 * (k - 1) + c] = blk[i][c]
        return L

    def K2(t):
        return [[sum((vm.Kel(Fr(re), Fr(im)) * t ** int(k) for k, re, im in data['K2'][i][j]), vm.ZERO) * vm.Kel(1 / (1 + t * t))
                 for j in range(12)] for i in range(16)]
    good = True
    degK = max(int(k) for row in data['K2'] for e in row for k, re, im in e if Fr(re) != 0 or Fr(im) != 0)
    for t in [Fr(0), Fr(1, 3), Fr(-2), Fr(5, 7), Fr(3)]:
        good &= vm.iszero_m(vm.mmul(Lsym([Fr(1), Fr(0), t, Fr(0)]), K2(t)))
    good &= degK <= 3
    nz = lambda e: any(Fr(re) != 0 or Fr(im) != 0 for k, re, im in e)
    one_plus_t2 = lambda e: sorted((int(k), Fr(re), Fr(im)) for k, re, im in e if Fr(re) != 0 or Fr(im) != 0) == [(0, Fr(1), Fr(0)), (2, Fr(1), Fr(0))]
    outside = [i for i in range(16) if i not in (4, 5, 6, 7)]
    pos = [[j for j in range(12) if nz(data['K2'][i][j])] for i in outside]
    idr = all(len(p_) == 1 and one_plus_t2(data['K2'][i][p_[0]]) for i, p_ in zip(outside, pos)) and sorted(p_[0] for p_ in pos) == list(range(12))
    rec(f'K\'(t): L(1,0,t,0) (1+t^2)K\'(t) == 0 (entries of (1+t^2)K\' have degree {degK} <= 3, L linear, so degree <= 4: checked at 5 values of t) '
        "and K' restricted to the rows outside block 2 is the 12 x 12 identity (full column rank for every t)", good and idr)
    rec('recorded M2 degrees <= (6, 4)', all(max((i for i, j in p), default=0) <= 6 and max((j for i, j in p), default=0) <= 4 for r in M2 for p in r))
    rec('recorded M1 degrees <= (4, 4)', all(max((i for i, j in p), default=0) <= 4 and max((j for i, j in p), default=0) <= 4 for r in M1 for p in r))
    K1c = data['K1']
    g1 = g2 = True
    for e in [Fr(-4000), Fr(-700, 3), Fr(0), Fr(1500, 7), Fr(9000)]:
        S = [vm.madd(Srat[r], vm.msc(vm.Kel(e / 35736), Sq[r])) for r in range(5)]
        M, Bj = vm.build_parts(G, g3, S); Minv = vm.inverse(M)
        P1, P3 = vm.mmul(Minv, Bj[1]), vm.mmul(Minv, Bj[3])
        for tau in [Fr(1, 50), Fr(1, 10), Fr(2, 7), Fr(1, 2), Fr(5, 6)]:                       # ray 1
            s = 2 * tau / (1 - tau * tau); x = (1 + tau * tau) / (1 - tau * tau)
            Dn = (tau * tau - 1) ** 2 * (4 * e * e - 27625 ** 2) * (16 * e * e - 63361 ** 2)
            Kt = [[vm.Kel(Fr(q[0]), Fr(q[1])) + vm.Kel(Fr(q[2]), Fr(q[3])) * s for q in row] for row in K1c]
            PK = vm.mmul(vm.madd(vm.madd(vm.msc(s, P1), P3), vm.msc(x, vm.eye(16)), -vm.ONE), Kt)
            for i in range(16):
                for j in range(12):
                    z = PK[i][j] * Dn; rv, iv = bt.poly_eval(M1[i][j], tau, e)
                    g1 &= z.a == rv and z.b == iv and z.c == 0 and z.d == 0
        for sg in [Fr(0), Fr(1, 9), Fr(1, 4), Fr(2, 5), Fr(1, 2), Fr(2, 3), Fr(4, 5)]:            # ray 2 (7 values)
            t = 2 * sg / (1 - sg * sg); x = (1 + sg * sg) / (1 - sg * sg)
            Dn = (sg * sg - 1) * (sg * sg + 1) ** 2 * (4 * e * e - 27625 ** 2) * (16 * e * e - 63361 ** 2)
            PK = vm.mmul(vm.madd(vm.madd(P1, vm.msc(t, P3)), vm.msc(x, vm.eye(16)), -vm.ONE), K2(t))
            for i in range(16):
                for j in range(12):
                    z = PK[i][j] * Dn; rv, iv = bt.poly_eval(M2[i][j], sg, e)
                    g2 &= z.a == rv and z.b == iv and z.c == 0 and z.d == 0
    rec('M1 == independently computed Dn (P - |xi|) K on the 5 x 5 grid', g1)
    rec('M2 == independently computed Dn\' (P - |xi\'|) K\' on the 7 x 5 grid', g2)
    boxes = data['boxes']; nr = nc = ne = 0
    for b in boxes:
        tl, th = Fr(b['t_lo']), Fr(b['t_hi']); el, eh = Fr(b['eta_lo']), Fr(b['eta_hi'])
        vl, vh = Fr(b['v_lo']), Fr(b['v_hi']); elo, ehi = Fr(b['eps_lo']), Fr(b['eps_hi'])
        C = [[(Fr(c[0]), Fr(c[1])) for c in row] for row in b['C']]
        Mp = M1 if b['ray'] == 1 else M2
        nr += vx.rank8_subdivided(Mp, b['rows'], b['cols'], C, tl, th, el, eh)
        nc += bt.segment_containment(tl, th, vl, vh, branch=b['branch'], ray=b['ray'], touch_zero=b['touch_zero'])
        a_, z_ = (el, eh) if b['sign'] > 0 else (-eh, -el)
        ne += (vl == elo * elo * D and vh == ehi * ehi * D and 0 <= a_ and a_ * a_ <= vl and z_ * z_ >= vh
               and abs(el) < Fr(27625, 2) and abs(eh) < Fr(27625, 2) and th < 1)
    rec(f'rank-8 Taylor test on all {len(boxes)} boxes ({nr} passed)', nr == len(boxes))
    rec(f'containment / branch identity on all boxes ({nc} passed)', nc == len(boxes))
    rec(f'eta enclosures and v = eps^2 D on all boxes ({ne} passed)', ne == len(boxes))
    eL, eU = Fr(data['end_eps']['lower']), Fr(data['end_eps']['upper'])

    def above(e, a, b, n):          # e^2 >= (a - b sqrt(n))/823  <=>  b sqrt(n) >= a - 823 e^2
        rhs = a - 823 * e * e
        return rhs <= 0 or b * b * n >= rhs * rhs
    rec('end eps of the lower end boxes >= eps_L', above(eL, 34404, 12, 5659689))
    rec('end eps of the upper end boxes >= eps_U', above(eU, 58228, 20, 1645345))
    for ray, br, e0, e1 in [(1, 'lower', Fr(266, 100), eL), (2, 'upper', Fr(628, 100), eU)]:
        for sg_ in (1, -1):
            seg = [b for b in boxes if b['ray'] == ray and b['sign'] == sg_]
            rec(f'ray {ray} {br} sign {sg_:+d}: contiguous {e0} .. {e1}, last box touches 0',
                Fr(seg[0]['eps_lo']) == e0 and Fr(seg[-1]['eps_hi']) == e1 and seg[-1]['touch_zero']
                and all(Fr(seg[i]['eps_hi']) == Fr(seg[i + 1]['eps_lo']) for i in range(len(seg) - 1)))
    line = 'SECOND CHECK OF THE ENDPOINT-GAP CERTIFICATE: ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if out:
        open(out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
