"""verify_m2_extended_fractions.py -- second check of the extended Milestone 2 certificate (pure fractions; no SymPy).

Usage:  python3 -I verify_m2_extended_fractions.py data/m2_extended_certificate.json [--out LOG]

  1. Matrix identity. M(tau, eta) = Dn (P(xi) - |xi|) K is recomputed from this project's independent fractions
     construction of P (verify_m2_fractions.model / build_parts, no SymPy, no rs_model) at a 5 x 5 grid of rational
     (tau, eta) and compared entry by entry with the recorded polynomials. Each recorded entry has degree <= 4 in tau
     and in eta (checked), and so does the true entry (S^rho affine in eta, det S^0 of degree 4 cleared by Dn; P, K
     affine in s, (1 - tau^2)^2 in Dn), so agreement on the grid proves the recorded matrix.
  2. Rank-8 test on every box: the same mathematical bound as the planner (Taylor expansion at the box centre,
     |p - p(c)| <= sum |c'_ij| ha^i hb^j) but an independent implementation: the Taylor shift is computed by repeated
     synthetic division (Horner), not by binomial expansion as in m2_box_tools, and the box may be subdivided (up to
     4 levels) with the recorded C kept. (A first version used naive interval arithmetic on monomials; it could not
     confirm 46 of 108 boxes because of the wide eta ranges, which is a weakness of that bound, not a failed claim;
     see notes/FAILED_ATTEMPTS.md F-6.)
  3. Containment and branch identity (m2_box_tools.segment_containment, exact Sturm counts), eta enclosures
     (eta range covers [sqrt(v_lo), sqrt(v_hi)] with the right sign), v = eps^2 * 2450894, contiguity of segments from
     |eps| = 1 to the targets, and targets below the exact thresholds eps_L, eps_U.
"""
import hashlib
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_m2_fractions as vm  # noqa: E402
import m2_box_tools as bt  # noqa: E402

D = 2450894


def ipow(lo, hi, n):
    """Interval [lo, hi]^n."""
    if n == 0:
        return Fr(1), Fr(1)
    cands = [lo ** n, hi ** n]
    if lo < 0 < hi and n % 2 == 0:
        return Fr(0), max(cands)
    return min(cands), max(cands)


def imul(a, b):
    ps = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return min(ps), max(ps)


def enclose(poly, tl, th, el, eh):
    """Enclosures (re_lo, re_hi), (im_lo, im_hi) of poly over the box."""
    rl = rh = il = ih = Fr(0)
    for (i, j), (cr, ci) in poly.items():
        m = imul(ipow(tl, th, i), ipow(el, eh, j))
        for c, which in ((cr, 'r'), (ci, 'i')):
            if c == 0:
                continue
            lo, hi = (c * m[0], c * m[1]) if c > 0 else (c * m[1], c * m[0])
            if which == 'r':
                rl += lo; rh += hi
            else:
                il += lo; ih += hi
    return (rl, rh), (il, ih)


def rank8_interval(Mpoly, rows, cols, C, tl, th, el, eh):
    enc = {(r, c): enclose(Mpoly[r][c], tl, th, el, eh) for r in rows for c in cols}
    for i in range(8):
        rs = Fr(0)
        for j in range(8):
            re_lo = re_hi = im_lo = im_hi = Fr(0)
            for k in range(8):
                (ar, br_), (ai, bi) = enc[(rows[k], cols[j])]
                cr, ci = C[i][k]
                # (cr + i ci)(x + i y), x in [ar, br_], y in [ai, bi]
                for coef, iv, tgt in ((cr, (ar, br_), 'rr'), (-ci, (ai, bi), 'rr'), (cr, (ai, bi), 'ii'), (ci, (ar, br_), 'ii')):
                    lo, hi = (coef * iv[0], coef * iv[1]) if coef >= 0 else (coef * iv[1], coef * iv[0])
                    if tgt == 'rr':
                        re_lo += lo; re_hi += hi
                    else:
                        im_lo += lo; im_hi += hi
            if i == j:
                re_lo -= 1; re_hi -= 1
            rs += max(abs(re_lo), abs(re_hi)) + max(abs(im_lo), abs(im_hi))
            if rs >= 1:
                return False
    return True


def taylor_shift_1d(coeffs, c):
    """Coefficients of q(c + a) from q(x) (low degree first), by repeated synthetic division."""
    a = list(coeffs); n = len(a); out = []
    for k in range(n):
        # divide by (x - c): remainder is q^(k)(c)/k!
        rem = Fr(0); new = [Fr(0)] * (len(a) - 1) if len(a) > 1 else []
        acc = Fr(0)
        for i in range(len(a) - 1, -1, -1):
            acc = acc * c + a[i]
            if i > 0:
                new[i - 1] = acc
        out.append(acc)
        a = new
        if not a:
            break
    out += [Fr(0)] * (n - len(out))
    return out


def shift2(poly, tc, ec):
    """{(i, j): (re, im)} of p(tc + a, ec + b), via 1-D shifts in eta then tau."""
    di = max(i for i, j in poly); dj = max(j for i, j in poly)
    res = {}
    for part in (0, 1):
        grid = [[poly.get((i, j), (Fr(0), Fr(0)))[part] for j in range(dj + 1)] for i in range(di + 1)]
        grid = [taylor_shift_1d(row, ec) for row in grid]                      # shift in eta
        for j in range(dj + 1):
            col = taylor_shift_1d([grid[i][j] for i in range(di + 1)], tc)      # shift in tau
            for i in range(di + 1):
                if col[i] != 0:
                    o = res.get((i, j), [Fr(0), Fr(0)]); o[part] = col[i]; res[(i, j)] = o
    return {k: (v[0], v[1]) for k, v in res.items()}


def rank8_taylor(Mpoly, rows, cols, C, tl, th, el, eh):
    tc, ec, ha, hb = (tl + th) / 2, (el + eh) / 2, (th - tl) / 2, (eh - el) / 2
    sh = {(r, c): shift2(Mpoly[r][c], tc, ec) for r in rows for c in cols}
    for i in range(8):
        rs = Fr(0)
        for j in range(8):
            ent = {}
            for k in range(8):
                cr, ci = C[i][k]
                for m, (pr, pi) in sh[(rows[k], cols[j])].items():
                    o = ent.get(m, (Fr(0), Fr(0)))
                    ent[m] = (o[0] + cr * pr - ci * pi, o[1] + cr * pi + ci * pr)
            c0 = ent.get((0, 0), (Fr(0), Fr(0)))
            rs += abs(c0[0] - (1 if i == j else 0)) + abs(c0[1])
            rs += sum((abs(v[0]) + abs(v[1])) * ha ** a * hb ** b for (a, b), v in ent.items() if (a, b) != (0, 0))
            if rs >= 1:
                return False
    return True


def rank8_subdivided(Mpoly, rows, cols, C, tl, th, el, eh, depth=4):
    if rank8_taylor(Mpoly, rows, cols, C, tl, th, el, eh):
        return True
    if depth == 0:
        return False
    tm, em = (tl + th) / 2, (el + eh) / 2
    return all(rank8_subdivided(Mpoly, rows, cols, C, a, b, c, d, depth - 1)
               for a, b in ((tl, tm), (tm, th)) for c, d in ((el, em), (em, eh)))


def main():
    path = sys.argv[1]; out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    raw = open(path, 'rb').read(); data = json.loads(raw); log = []; ok_all = True

    def say(s=''):
        print(s, flush=True); log.append(s)

    def rec(name, good):
        nonlocal ok_all
        ok_all &= bool(good); say(f'   {"OK  " if good else "FAIL"} {name}: {bool(good)}')
    import platform
    say('Second check of the extended Milestone 2 certificate (fractions only, no SymPy)')
    say(f'python {platform.python_version()} | input {os.path.basename(path)} sha256 {hashlib.sha256(raw).hexdigest()}')
    say(f'script sha256 {hashlib.sha256(open(__file__, "rb").read()).hexdigest()}')
    say(f'm2_box_tools.py sha256 {hashlib.sha256(open(os.path.join(HERE, "m2_box_tools.py"), "rb").read()).hexdigest()}')
    say(f'verify_m2_fractions.py sha256 {hashlib.sha256(open(os.path.join(HERE, "verify_m2_fractions.py"), "rb").read()).hexdigest()}')
    Mpoly = [[bt.poly_from_json(e) for e in row] for row in data['M']]
    Kc = data['K']
    say('1. Matrix identity on a 5 x 5 grid')
    rec('every recorded entry has degree <= 4 in tau and in eta',
        all(max((i for i, j in p), default=0) <= 4 and max((j for i, j in p), default=0) <= 4 for row in Mpoly for p in row))
    G, g3, Srat, Sq, cancel = vm.model()
    rec('S_mu - (16/3) i S_q == 0', cancel)
    good = True
    for e in [Fr(-3000), Fr(-700, 3), Fr(0), Fr(1500, 7), Fr(9000)]:
        S = [vm.madd(Srat[r], vm.msc(vm.Kel(e / 35736), Sq[r])) for r in range(5)]
        M, Bj = vm.build_parts(G, g3, S); Minv = vm.inverse(M)
        for t in [Fr(1, 50), Fr(1, 10), Fr(2, 7), Fr(1, 2), Fr(5, 6)]:
            s = 2 * t / (1 - t * t); x = (1 + t * t) / (1 - t * t)
            Dn = (t * t - 1) ** 2 * (4 * e * e - 27625 ** 2) * (16 * e * e - 63361 ** 2)
            Kt = [[vm.Kel(Fr(q[0]), Fr(q[1])) + vm.Kel(Fr(q[2]), Fr(q[3])) * s for q in row] for row in Kc]
            NB = vm.madd(vm.madd(vm.msc(s, Bj[1]), Bj[3]), vm.msc(x, M), -vm.ONE)
            PK = vm.mmul(Minv, vm.mmul(NB, Kt))
            for i in range(16):
                for j in range(12):
                    rv, iv = bt.poly_eval(Mpoly[i][j], t, e)
                    z = PK[i][j] * Dn
                    if not (z.a == rv and z.b == iv and z.c == 0 and z.d == 0):
                        good = False
    rec('recorded M(tau, eta) == independently computed Dn (P - |xi|) K at all 25 grid points, all 192 entries', good)
    say('2-3. Boxes')
    boxes = data['boxes']; targets = {k: Fr(v) for k, v in data['targets'].items()}
    nrank = ncont = nenc = 0
    for b in boxes:
        tl, th = Fr(b['t_lo']), Fr(b['t_hi']); el, eh = Fr(b['eta_lo']), Fr(b['eta_hi'])
        vl, vh = Fr(b['v_lo']), Fr(b['v_hi']); elo, ehi = Fr(b['eps_lo']), Fr(b['eps_hi'])
        C = [[(Fr(c[0]), Fr(c[1])) for c in row] for row in b['C']]
        nrank += rank8_subdivided(Mpoly, b['rows'], b['cols'], C, tl, th, el, eh)
        ncont += bt.segment_containment(tl, th, vl, vh, branch=b['branch'])
        a_, z_ = (el, eh) if b['sign'] > 0 else (-eh, -el)
        nenc += (vl == elo * elo * D and vh == ehi * ehi * D and 0 <= a_ and a_ * a_ <= vl and z_ * z_ >= vh
                 and 1 <= elo < ehi and abs(el) < Fr(27625, 2) and abs(eh) < Fr(27625, 2))
    rec(f'rank-8 Taylor test (independent shift, subdivided if needed) passes on all {len(boxes)} boxes ({nrank} passed)', nrank == len(boxes))
    rec(f'containment + branch identity on all boxes ({ncont} passed)', ncont == len(boxes))
    rec(f'eta enclosures, v = eps^2 * 2450894, Dn != 0 on all boxes ({nenc} passed)', nenc == len(boxes))
    for name in ('lower', 'upper'):
        for sg in (1, -1):
            seg = [b for b in boxes if b['branch'] == name and b['sign'] == sg]
            cov = bool(seg) and Fr(seg[0]['eps_lo']) == 1 and Fr(seg[-1]['eps_hi']) == targets[name] and all(
                Fr(seg[i]['eps_hi']) == Fr(seg[i + 1]['eps_lo']) for i in range(len(seg) - 1))
            rec(f'{name}, sign {sg:+d}: contiguous cover of [1, {targets[name]}]', cov)
    # thresholds: tl^2 < (34404 - 12 sqrt(5659689))/823  <=>  12 sqrt(5659689) < 34404 - 823 tl^2 (rhs > 0) <=> squares
    def below(t, a, b, n):          # t^2 < (a - b sqrt(n))/823
        rhs = a - 823 * t * t
        return rhs > 0 and b * b * n < rhs * rhs
    rec('lower target^2 < eps_L^2 = (34404 - 12 sqrt(5659689))/823', below(targets['lower'], 34404, 12, 5659689))
    rec('upper target^2 < eps_U^2 = (58228 - 20 sqrt(1645345))/823', below(targets['upper'], 58228, 20, 1645345))
    line = 'SECOND CHECK OF EXTENDED MILESTONE 2 CERTIFICATE: ' + ('ALL PASS' if ok_all else 'SOME CHECK FAILED')
    say(line)
    if out:
        open(out, 'w', encoding='utf-8').write('\n'.join(log) + '\n')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
