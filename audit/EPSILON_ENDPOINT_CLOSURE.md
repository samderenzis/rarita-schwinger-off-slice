# EPSILON_ENDPOINT_CLOSURE — closing 2.66 ≤ |ε| ≤ ε_L and 6.28 ≤ |ε| ≤ ε_U, endpoints included

Date: 2026-11-22 (GEN43 Phase 1). Author: Claude (AI). Review status: AI-checked only. Tags as in `OFFSLICE_COERCIVITY_FINAL.md`.

## 1. The gap
Before this step, rank M ≥ 8 (geometric multiplicity 4) was certified for the lower crossing up to |ε| = 2.66 and for the upper crossing up to |ε| = 6.28. The crossings exist up to ε_L ≈ 2.6674 and ε_U ≈ 6.2912 (exact values in JORDAN_STRUCTURE_COMPLETE). Two problems blocked a direct extension:
- near ε_L the lower crossing runs into s = 0, where the ordinary box containment (box strictly inside (0, 1) in τ) cannot be used;
- near ε_U the upper crossing runs to s = ∞, i.e. τ → 1, where the ray-1 matrix M has the factor (τ² − 1)² and degenerates.

## 2. Coverage after this step (both signs of ε)
| Branch | |ε| range | Ray | Certificate | Boxes | Worst ‖CA − I‖_∞ |
|---|---|---|---|---|---|
| both | (0, 1] | 1 | `data/m2_minor.json` (RESULT_m2_q_direction) | 1 minor | — (Taylor bound) |
| lower | [1, 2.66] | 1 | `data/m2_extended_certificate.json` | 20 | ≤ 0.9708 |
| upper | [1, 6.28] | 1 | `data/m2_extended_certificate.json` | 88 | 0.9708 |
| lower | [2.66, 1067/400] ∋ ε_L | 1 | `data/m2_gap_certificate.json` | 10 | 0.0866 |
| upper | [6.28, 62913/10000] ∋ ε_U | 2 | `data/m2_gap_certificate.json` | 12 | 0.1635 |

1067/400 = 2.6675 ≥ ε_L and 62913/10000 = 6.2913 ≥ ε_U are checked exactly (squares compared with the exact surds). The extended certificate was re-verified with the modified `m2_box_tools.py` (log `results/regression_verify_m2_extended_new_boxtools_2026-11-22.txt`, ALL PASS), so the box-tool change did not break the earlier certificate.

## 3. Method
**Ray 2 for the upper branch.** Use ξ′ = (1, 0, t, 0), t = 2σ/(1 − σ²), σ ∈ (0, 1). ξ′ is a positive multiple of (s, 0, 1, 0) with s = 1/t, so the Jordan structure at |ξ′| is the same as at |ξ| (P is linear in ξ). The upper crossing s → ∞ becomes t → 0, a regular point. In the crossing polynomial the roles of the constant and leading coefficients (c₀ and c₂ in u) swap; the upper branch is the smallest root σ in (0, 1).

**Basis K′ of ker L(ξ′).** The SymPy nullspace of L(ξ′) contains 1/t. Instead, K′ is built by hand: the identity on the 12 components outside block 2, with block 2 solved through B(t) = Γ⁰(Γ¹ + tΓ³)Γ², det B = (1 + t²)² ≠ 0 for every real t. So K′ has full column rank for all real t including t = 0, and (1 + t²)K′ is polynomial of degree 2. Second check: L(ξ′)(1 + t²)K′ = 0 at 5 values of t (degree ≤ 4 in t, so this is a proof), and the rows outside block 2 are exactly (1 + t²)·I₁₂.

**Ordinary boxes.** As in the extended certificate: Sturm containment of the branch in the box for the whole ε-segment (one root at the start, none on the box sides, branch identity), and ‖CA − I‖_∞ < 1 for an 8×8 submatrix A of M with a rational C, by exact Taylor bounds.

**End boxes (touch_zero mode).** The last box on each branch is [0, t_hi] × [ε_lo, ε_end] with ε_end beyond the threshold. Exact checks: at ε_lo, exactly one root of the crossing polynomial in (0, t_hi), not at 0 or t_hi; for all ε in the segment, no root at t = t_hi. The rank bound holds on the whole closed box, including t = 0.

*Why this suffices (continuity argument).* The u-discriminant of C is ∝ ε²(91809 − 823ε²), nonzero for 0 < ε² < 111.5, so the two crossing roots are simple and depend continuously on ε. On the end segment the branch starts inside (0, t_hi). It cannot cross t_hi (no root there). The only way out is through t = 0, which happens exactly at the threshold, because the constant coefficient (ray 1) or leading coefficient (ray 2) of C in u vanishes, within the segment, only at the threshold (its other zero is at ε ≈ 8.746 for ray 1 and ε ≈ 10.09 for ray 2). So for every ε in [ε_lo, threshold) the crossing lies in the box; at the threshold it is at t = 0, also in the box; past the threshold it does not exist and nothing is claimed. Roots that might enter the box from the complex plane do not matter: the rank bound holds on the whole box. **ESTABLISHED.**

## 4. Algebraic multiplicity at the endpoints
At ε = ε_L the crossing is ξ = (0, 0, 1, 0); at ε = ε_U it is ξ = (1, 0, 0, 0). The ray dispersion relation is written in s, so these need a separate check. `audit/endpoint_multiplicity.py`:
- Q has total degree ≤ 4 in (x, s); its homogenisation Q_h(x, s, b) satisfies Q_h(x, s, 1) = Q, and lc_x is unchanged;
- lower endpoint (s = 0, ε = ε_L): X = 1 is a root of Q with dQ/dX = −67603406660674704 − 27990756539088√5659689 ≈ −1.342·10¹⁷ ≠ 0;
- upper endpoint (b = 0, ε = ε_U): X = 1 is a root with derivative 422979956239154544 + 340756598051280√1645345 ≈ 8.601·10¹⁷ ≠ 0;
- the quotient factor does not vanish at r = 1 (1 − 1/4 ≠ 0).

So r = 1 has algebraic multiplicity exactly 5 on ker L at both endpoints. With rank M = 8 from the end boxes: J₂ ⊕ 3J₁ at the endpoints too. **ESTABLISHED (SymPy only).**

## 5. Margins
| Certificate | Boxes | row-sum min / median / max | Margin to 1 | τ-widths | ε-widths |
|---|---|---|---|---|---|
| extended | 108 | 0.087 / 0.636 / 0.971 | 0.029 (worst: upper, 223/40 ≤ ε ≤ 57/10) | 4.2·10⁻³ – 5.5·10⁻² | 2.9·10⁻³ – 0.5 |
| gaps | 22 | 0.011 / 0.038 / 0.164 | 0.836 | 4.2·10⁻⁴ – 2.4·10⁻² | 1.0·10⁻⁴ – 4.0·10⁻³ |

The extended certificate's worst margin (0.029) is small but exact; it is a strict rational inequality, not a floating-point estimate.

## 6. Second check and shared code
`scripts/verify_m2_gaps_fractions.py` (no SymPy): recomputes M1 on a 5×5 grid and M2 on a 7×5 grid from the independent fractions construction of P (so every entry is proved, given the degree bounds 4/4 and 6/4); checks K′; redoes every rank-8 bound with an independent Taylor-shift implementation; re-checks containment, η enclosures, end ε ≥ thresholds and contiguity. ALL PASS. The containment routine (`m2_box_tools.segment_containment`) is shared with the planner.

A bug in a first version of this checker counted zero polynomial entries as nonzero (FAILED_ATTEMPTS F-7); it was fixed and the degree-4 test was corrected to use 5 values of t.

## 7. Status
Gap closure: **ESTABLISHED**, with a second check that shares the containment code. Endpoint multiplicity: **ESTABLISHED (SymPy only)**. No interval claim is made beyond the thresholds; for |ε| > ε_L (lower) and > ε_U (upper) the crossings are not real covectors.
