# JORDAN_STRUCTURE_COMPLETE — exact Jordan classification on the full crossing intervals

Date: 2026-11-22 (GEN43 Phase 1 final form; supersedes `JORDAN_STRUCTURE_CLASSIFICATION.md`, which stated the result on 0 < |ε| ≤ 2.66 / 6.28 and is kept unchanged). Author: Claude (AI). Review status: AI-checked only. Tags as in `OFFSLICE_COERCIVITY_FINAL.md`.

## Setting
Witness m² = 823/142944, λ² = 1/2, b₁ = 1, b₂ = 1/3, any π′, q = (1 + ε)q_slice. Exact thresholds:
ε_L = √((34404 − 12√5659689)/823) ≈ 2.6674 and ε_U = √((58228 − 20√1645345)/823) ≈ 6.2912.
For 0 < |ε| < ε_L the lower crossing exists at a real ratio u₋(ε) ∈ (0, ∞); at |ε| = ε_L it reaches u = 0, i.e. ξ ∝ (0,0,1,0). For 0 < |ε| < ε_U the upper crossing exists; at |ε| = ε_U it reaches u = ∞, i.e. ξ ∝ (1,0,0,0). "Crossing covector" below includes these limit directions at the endpoints.

## Theorem J (complete form)
Let ξ be a crossing covector for the lower branch with 0 < |ε| ≤ ε_L, or for the upper branch with 0 < |ε| ≤ ε_U. Then r = |ξ| is an eigenvalue of A = P(ξ)|ker L(ξ) with algebraic multiplicity exactly 5 and geometric multiplicity exactly 4. Its Jordan form on the generalised eigenspace is

  **J₂(r) ⊕ J₁(r) ⊕ J₁(r) ⊕ J₁(r).**

The same holds at −r (P(−ξ) = −P(ξ)). A Jordan chain v → w in ker L exists, and v spans the one-dimensional space (A − r) ker(A − r)².

*Proof.*
1. alg ≥ 5: RESTRICTED_SYMBOL_COMPLETE_PROOF R6–R7 (crossing ⇒ C = 0).
2. alg ≤ 5: interior points by Proposition RI (RESULTANT_AND_ROOT_ISOLATION §6). Endpoints by the homogenised multiplicity check (EPSILON_ENDPOINT_CLOSURE §4).
3. geo ≥ 4 everywhere: rank M ≤ 8 identically (R8).
4. geo ≤ 4, i.e. rank M ≥ 8, at every crossing in the stated ranges:
   - 0 < |ε| ≤ 1, both branches, ray 1: one 8×8 minor with a Taylor bound (`RESULT_m2_q_direction`, `data/m2_minor.json`);
   - 1 ≤ |ε| ≤ 2.66 (lower) and 1 ≤ |ε| ≤ 6.28 (upper), ray 1: 108 boxes (`data/m2_extended_certificate.json`, worst bound 0.9708);
   - 2.66 ≤ |ε| ≤ 1067/400 (≥ ε_L), lower, ray 1, and 6.28 ≤ |ε| ≤ 62913/10000 (≥ ε_U), upper, ray 2: 22 boxes (`data/m2_gap_certificate.json`, worst bound 0.1635); the end boxes contain the endpoint directions s = 0 and t = 0.
   In each box the crossing is shown to lie inside the box by exact Sturm containment with branch identification.
5. Rotation equivariance moves the ray results to every covector on the crossing cones (SymPy only).
6. Blocks: number = geo = 4, sizes sum to alg = 5 ⇒ 2, 1, 1, 1. ∎

**Status: ESTABLISHED.** Weakest links: the resultant step, the endpoint multiplicity and rotation equivariance are SymPy only; the box containment routine is shared between the box planner and the second checker.

## What changes where
| Locus | Structure at the light eigenvalue on ker L | Status |
|---|---|---|
| ε = 0 (slice) | alg 6, geo 4: J₂ ⊕ J₂ ⊕ J₁ ⊕ J₁ (both crossings merge) | ESTABLISHED (parent) |
| 0 < \|ε\| ≤ ε_L (lower), 0 < \|ε\| ≤ ε_U (upper) | J₂ ⊕ 3J₁ at each crossing | ESTABLISHED (Theorem J) |
| off the crossing cones | semisimple, alg = geo = 4 | ESTABLISHED |
| ε_L < \|ε\| (lower), ε_U < \|ε\| < ≈8.746 | no real light/G crossing (all three coefficients of C in u negative) | EXPLORATORY |
| ≈8.746 < \|ε\| < 8.823 | lower crossing reappears | EXPLORATORY |
| ε ≈ 8.823, 10.118 | det S⁰ = 0: elimination of ψ₀ fails; outside the framework | ESTABLISHED (exact roots) |
| ε ≈ 10.133, 10.397 | real triple crossings, alg = 6, type unknown | EXPLORATORY |

## Answers to the Phase 1 questions
- Is there exactly one J₂? **Yes**, at every crossing in the full existence intervals near the slice.
- Can J₃ occur there? **No** (alg = 5, geo = 4 forces sizes 2,1,1,1).
- Can two J₂ blocks occur there? **No**; two J₂ occur only at ε = 0.
- Where could the structure change? Only where disc_X Q and C share a real root with u > 0 (the two triple crossings above) or where det S⁰ = 0, both far outside the intervals.

## Strongest justified conclusion
For the witness parameters, the slice's J₂ ⊕ J₂ splits for every ε ≠ 0 into one J₂ at each of two separate crossing cones, and this persists exactly, with no further degeneration, until each crossing leaves the real covectors at ε_L and ε_U. This is a statement about the constrained principal symbol on a one-parameter family. It is not an open-set theorem in (m, λ, b, q), and it carries no physical, causality or hyperbolicity claim.
