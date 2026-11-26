# Plan: is the witness isolated?

This plan follows the parent project's [ROADMAP §2](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/main/ROADMAP.md#2-is-the-witness-isolated-or-part-of-a-family). If the two ever disagree, the ROADMAP is the reference until this file is updated and the change is logged in CHANGELOG.md.

## Background that must not be assumed
- On the Pauli slice, the parent records report defect 2 at every crossing on 27 tested slices.
- Off the slice, the degenerate double eigenvalue splits at the tested values. Whether a defect survives is an **open** question (parent erratum G-66).
- The on-slice factorisation and the branch names (light branch, G-branch) must not be assumed to continue off the slice.

## Sequence
1. **Make q independent.** The parent scripts `witness_min.py` and `bhmap_defs.py` hard-wire q = λ√(2/3)·m. The factor `kap` scales λ in both q and the Pauli term, so it does not leave the slice.
2. **Re-verify the reduction with q free.** Check that the divergence identity (paper, Prop. 3.2), the constraint intertwining, and the π′-factorisation P_{π′} − P_{π″} = (π′ − π″)M⁻¹TL hold, wherever det S⁰ ≠ 0 and M is invertible.
3. **Parametrise the family**, for example q = (1 + ε)·λ√(2/3)·m.
4. **Exploratory step:** for a few small values of ε, in high precision, factor the characteristic polynomial of P(ξ) and locate real crossings near the witness.
5. **Exact step:** at the chosen points, give the crossing covector as an exact algebraic number. Compute the ranks of (P − r)^j on ker L over that number field. If a defect is found, write down an explicit chain v → w inside ker L.

## Milestones
- **Milestone 1:** one exact off-slice constrained Jordan chain. Every hypothesis of the analytic argument is re-checked at that point:
  - det S⁰ ≠ 0;
  - the reduction identities;
  - coercivity on some ball;
  - the curved-symbol identification.

  This would be a second pointwise witness, not a family theorem.
- **Milestone 2:** a parametrised argument showing where such chains exist, and where the hypotheses hold (with uniformity, if a common ball is wanted). Only this would support a family-level statement.

## Progress (updated 2026-11-25; not yet released; AI-checked only)
- Steps 1–5 done. **Milestone 1 reached**, pointwise at four exact points: `results/RESULT_offslice_points_2026-11-22.md`.
- **Milestone 2 reached in one direction only:** for every 0 < |ε| ≤ 1 (q varies; m, λ, b₁, b₂ fixed at the witness values), both crossings carry a Jordan chain on ker L, and the hypotheses hold at each ε: `results/RESULT_m2_q_direction_2026-11-22.md`. Extended to the lower crossing for 0 < |ε| ≤ 2.66 and the upper for 0 < |ε| ≤ 6.28, almost their whole existence ranges (ε_L ≈ 2.6674, ε_U ≈ 6.2912): `results/RESULT_m2_extended_2026-11-22.md`. The gaps are now closed up to and including ε_L and ε_U, and the Jordan form is exactly J₂ ⊕ 3J₁ there (`audit/JORDAN_STRUCTURE_COMPLETE.md`, `audit/EPSILON_ENDPOINT_CLOSURE.md`).
- Beyond the family (2026-11-25): at the four real triple crossings (ε ≈ ±10.133, ±10.397), the Jordan type is J₃ ⊕ 3J₁, certified externally by the GEN43 engine: `results/RESULT_triple_crossings_GEN43_2026-11-23.md`. The owner re-ran every script in this repository on his own machine: ALL MATCH (`results/OWNER_REPRODUCTION_2026-11-25.md`).
- Still open:
  - varying m, λ, b₁, b₂ (an open set in parameter space);
  - uniformity in ε of the radius and constants (needed for a common ball);
  - (done: written coercivity proof with explicit constants, `audit/OFFSLICE_COERCIVITY_FINAL.md`)
  - independent (non-AI) review.

## Outcomes, each reported separately
- **Positive:** an exact chain at a tested crossing, leading to Milestone 1.
- **Negative:** the tested crossings are exactly semisimple. This rules out the mechanism at those crossings only.
- **No crossing found:** this is not a proof that no crossing exists. That would need an exact exclusion argument.

## Not in scope
- Persistence on an open set, unless separately proved.
- Which of existence, uniqueness or the bound fails.
- Physics.
