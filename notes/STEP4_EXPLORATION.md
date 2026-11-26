# Step 4 exploration: where the off-slice crossings are (2026-11-22)

Script: `scripts/step4_dispersion_and_crossings.py`; log: `results/step4_dispersion_and_crossings_2026-11-22.txt`.
Parts A–C and E are exact computations; part D is 60-digit numerics (exploration only). Witness m, λ, b₁, b₂; π′ = 1/2; q = (1 + ε)·q_slice.

## What was found
1. **q enters P only through η = ε·√2450894 = 35736·(q − q_slice).** P has entries in ℚ(i)(η). The irrational √2450894 cancels, as it does on the slice.
2. **Exact dispersion relation on ker L along ξ = (s, 0, 1, 0), symbolic ε:**
   charpoly(P|ker L) ∝ (x² − s² − 1)⁴ · Q(x², s², ε²). Q is quadratic in x², even in ε, with rational coefficients; at ε = 0 it is the square of the parent G-sheet factor. So the light branch x² = |ξ|² keeps multiplicity 4 off the slice, and the doubled G-sheet splits into two sheets.
3. **Crossings.** The light branch meets a Q-sheet where C(u, ε) = 0, u = s². C is quadratic in u. At ε = 0 it has the double root u* = 400/9801 (the parent crossing). Its discriminant in u is 95227628033791228445188816896·ε²·(91809 − 823ε²), positive for 0 < |ε| < 10.56. So **off the slice the parent crossing splits into two distinct real crossings**, moving apart roughly in proportion to ε.
4. **Numerics (exploration only).** At ε = 1/10 and 1/100, at both crossings, the restriction A of P to ker L has numerical dim ker(A − r)^j = 4, 5, 5: one Jordan block each. At ε = 0: 4, 6, 6 (two blocks). So the parent's two blocks seem to separate onto the two crossings.
5. **Choosing an exact point.** In terms of ε², the crossing discriminant is c(u + 1)(9u + 1)²(82875u + 1267), with c = 2¹⁰3³823²1489⁵. Points with rational ε², rational covector and rational eigenvalue come from rational points on the conic w² = 4467(u + 1)(82875u + 1267) with u and 1 + u sums of two rational squares (search near u*). The smallest such point found is **point A**, on the lower crossing. **Point B** (third by size) is on the upper crossing. A− and B− are the same crossings with ε of the opposite sign. See `scripts/points.py`; the exact checks are in step 5 and in `results/RESULT_offslice_points_2026-11-22.md`.

## Used only to find the point, not to check it
- The search assumes that the characteristic polynomial depends on ξ only through |ξ_a|² and |ξ_b|² (rotation invariance on each sphere). Step 5 does not use this: it computes directly at the actual covector.
