# Jordan structure at the off-slice crossings: exact classification on the certified ranges

Date: 2026-11-22. Author: Claude (AI). Review status: AI-checked only. Tags as in `OFFSLICE_COERCIVITY_AUDIT.md`.

## Question
At the certified crossings, is there exactly one J₂ block, with all others of size 1? Can J₃, or two J₂, occur? Where does the structure change?

## Theorem J (classification on the certified ranges)
Witness m, λ, b; any π′; q = (1 + ε)q_slice. Let ξ = (s, 0, 1, 0) be a crossing (C(s², ε²) = 0, s > 0) with either
- the lower crossing and 0 < |ε| ≤ 2.66, or
- the upper crossing and 0 < |ε| ≤ 6.28.

Then r = |ξ| is an eigenvalue of A = P(ξ)|ker L(ξ) with **algebraic multiplicity exactly 5** and **geometric multiplicity exactly 4**. So the Jordan form of A on the generalised r-eigenspace is

  **J₂(r) ⊕ J₁(r) ⊕ J₁(r) ⊕ J₁(r).**

There is exactly one Jordan block of size 2, there is no block of size ≥ 3, and there are never two blocks of size 2. A chain v → w in ker L exists, and v is unique up to scaling: it spans (A − r)·ker(A − r)², which is one-dimensional. By rotation equivariance (SymPy only), the same holds at every covector on the corresponding crossing cones. By ξ → −ξ, it also holds at the eigenvalue −|ξ|.

*Proof.*
1. **alg ≥ 5** by RESTRICTED_SYMBOL_PROOF R6–R7.
2. **alg ≤ 5:** we need disc_X Q ≠ 0 at every such crossing. `audit/jordan_classification.py` (SymPy, exact) shows:
   - Res_u(disc_X Q, C) = const·v²·p₅(v)², with p₅ irreducible of degree 5. Any common zero (u, v) of disc and C is a zero of this resultant, because the resultant lies in the ideal (disc, C).
   - v = ε² > 0 rules out the factor v². By exact root isolation, p₅ has no root in (0, 2.66²], so the lower range is clear.
   - On the upper range p₅ has exactly one root v₀ ∈ [187607633/5809776, 11905019/368671] (ε₀ ≈ 5.6826). For v in that interval the upper crossing u₊(v) lies in [17/10, 19/10]. This is shown by exact Sturm counts: one root at the left end, none above, none on the sides. On [17/10, 19/10] × [v₀ interval], disc_X Q is enclosed in [3.4·10³⁶, 9.5·10³⁶], which excludes 0. The common zero at v₀ is the *lower* crossing continued to u ≈ −0.015 < 0, which is not a real covector.
3. **geo = 4:** rank (P − r)|ker L = rank M ≤ 8 identically (fractions engine, 37 × 37 grid) and ≥ 8 on the certified boxes (`RESULT_m2_q_direction`, `RESULT_m2_extended`).
4. **Block count:** the number of Jordan blocks for r equals geo = 4, and their sizes sum to alg = 5. So the sizes are 2, 1, 1, 1. ∎

**Status: ESTABLISHED.** The inputs R6 (dispersion) and R8 (rank ≤ 8, box bounds) are REPRODUCED or second-checked. The resultant and the v₀ exclusion are SymPy-exact only.

## Where the structure changes (exact facts, with exploratory context)
| Locus | What happens | Status |
|---|---|---|
| ε = 0 (slice) | disc_X Q = 0 (factor v): both G-sheets coincide; at the crossing alg = 6, geo = 4 → J₂ ⊕ J₂ ⊕ J₁ ⊕ J₁ (parent: dims 4, 6, 6) | ESTABLISHED (parent) |
| 0 < \|ε\| ≤ 2.66 (lower), 0 < \|ε\| ≤ 6.28 (upper) | J₂ ⊕ 3J₁ (Theorem J) | ESTABLISHED |
| 2.66 < \|ε\| < ε_L ≈ 2.6674 (lower), 6.28 < \|ε\| < ε_U ≈ 6.2912 (upper) | alg = 5 still holds on the lower gap (p₅ has no root below 32) and on the upper gap (no root of p₅ in (6.28², ε_U²)); geo = 4 is not certified there (rank ≥ 8 missing) | INCOMPLETE (alg part ESTABLISHED) |
| ε = ε_L | the lower crossing reaches s = 0 (ξ = (0,0,1,0)) and leaves the real covectors | ESTABLISHED (exact threshold) |
| ε = ε_U | the upper crossing reaches s = ∞ (ξ → (1,0,0,0) direction) | ESTABLISHED (exact threshold) |
| ε ≈ 5.6826 (root of p₅) | triple coincidence on the continued lower branch at u ≈ −0.015 (no real covector) | ESTABLISHED (exact), irrelevant to real covectors |
| ε ≈ 10.133, 10.397 (roots of p₅) | real triple crossings (u ≈ 11.59 and u ≈ 0.932): alg = 6; Jordan type unknown | EXPLORATORY (outside the family; det S⁰ ≠ 0 there, since its roots are at ≈ 8.823 and ≈ 10.118) |
| off the crossing cones | light eigenvalue semisimple (alg = geo = 4) | ESTABLISHED (RESTRICTED_SYMBOL_PROOF, consequence 1) |

## Remaining questions
- Jordan type at the real triple crossings (ε ≈ 10.13, 10.40), and beyond ε_U generally.
- An analytic (in ε) construction of the chain vectors v(ε), w(ε). Not attempted.
- Possible defects between the two G-sheets away from the light cone (where disc_X Q = 0 at real covectors, other than at crossings). Not yet classified.

## Strongest justified conclusion
On the certified ranges, the constrained symbol has exactly one Jordan block of size 2 at each crossing eigenvalue, with all other blocks at that eigenvalue of size 1. The parent's two J₂ blocks on the slice separate into one J₂ at each of the two crossings for ε ≠ 0. This is now a theorem on the certified ranges, not only a pattern seen at sample points.
