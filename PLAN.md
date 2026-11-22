# Plan: is the witness isolated?

This plan follows the parent project's [ROADMAP §2](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/main/ROADMAP.md#2-is-the-witness-isolated). If the two ever disagree, the ROADMAP is the reference until this file is updated and the change is logged in CHANGELOG.md.

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

## Outcomes, each reported separately
- **Positive:** an exact chain at a tested crossing, leading to Milestone 1.
- **Negative:** the tested crossings are exactly semisimple. This rules out the mechanism at those crossings only.
- **No crossing found:** this is not a proof that no crossing exists. That would need an exact exclusion argument.

## Not in scope
- Persistence on an open set, unless separately proved.
- Which of existence, uniqueness or the bound fails.
- Physics.
