# Result: Jordan chains along almost the whole real crossing set in the q-direction   (date: 2026-11-22)

This extends `RESULT_m2_q_direction_2026-11-22.md`, which covered 0 < |ε| ≤ 1. The setting is the same: q varies, and m, λ, b₁, b₂ stay at the witness values.

```
Outcome type: POSITIVE (constrained Jordan chain found), on a one-parameter family

Parameter family
  λ² = 1/2, m² = 823/142944, b₁ = 1, b₂ = 1/3, any π′;  q = (1 + ε)·q_slice.

Statement (algebraic part)
  Lower crossing ξ₋(ε) = (s₋(ε), 0, 1, 0): exists (s₋ real, > 0) exactly for 0 < ε² < ε_L²,
      ε_L² = (34404 − 12√5659689)/823, ε_L ≈ 2.6674 (at ε_L the constant coefficient of C(u, ε²) vanishes: s₋ → 0).
      Jordan chain in ker L certified for every 0 < |ε| ≤ 2.66.
  Upper crossing ξ₊(ε): exists exactly for 0 < ε² < ε_U², ε_U² = (58228 − 20√1645345)/823, ε_U ≈ 6.2912
      (the leading coefficient of C vanishes there: s₊ → ∞).
      Jordan chain in ker L certified for every 0 < |ε| ≤ 6.28.
  Both signs of ε. The same holds at every covector related to these by rotations within each sphere factor and by
  scaling. At each, r = |ξ| has algebraic multiplicity ≥ 5 and geometric multiplicity ≤ 4 on ker L.
  Uncertified gaps: 2.66 < |ε| < ε_L (lower) and 6.28 < |ε| < ε_U (upper). The step 4-style numerics show
  geometric multiplicity 4 there too, but that is exploration only.

Proof structure
  |ε| ≤ 1: RESULT_m2_q_direction (one 8×8 minor, Taylor box). 1 ≤ |ε| ≤ target: 108 segments (lower: 10 per sign,
  upper: 44 per sign), each with
    - a box [t_lo, t_hi] × [η range] in (τ, η), s = 2τ/(1 − τ²), η = ε√2450894;
    - containment of the branch for the whole segment: exact Sturm counts in v = η² = 2450894ε² (rational). The
      crossing polynomial has exactly one root τ in the box at the segment start and none on the box sides; branch
      identity: no root in (0, t_lo] for the lower branch, none in [t_hi, 1) for the upper;
    - rank ≥ 8 of M = Dn·(P − |ξ|)K on the box: for an 8×8 submatrix A and a rational C, ‖CA − I‖_∞ < 1 on the
      box (exact Taylor bounds of the degree-≤4 entries; worst row sum 0.971). So A is invertible on the box.
  Exact thresholds and targets: ε_L, ε_U from the coefficients of C; 2.66² < ε_L² and 6.28² < ε_U² exactly.
  Hypotheses for every such ε: det S⁰ ≠ 0 (|ε| < 8.82); reduction identities (step 2); curved-symbol identification
  for all covectors and charges (m2_curved_symbol_general.py); coercivity by proof reading. Each ε is a separate
  pointwise theorem; no uniformity.

Engines and cross-checks
  Planner: scripts/m2_family_extended.py (SymPy builds M exactly; boxes are checked with scripts/m2_box_tools.py in
  exact fractions; mpmath only to choose boxes, rows and columns and an approximate C).
  Second check: scripts/verify_m2_extended_fractions.py (no SymPy):
    - M recomputed from the independent fractions construction of P on a 5×5 grid, so all 192 entries are proved;
    - the rank test re-done with an independent Taylor-shift implementation (synthetic division; boxes may be
      subdivided);
    - containment, branch identity, η enclosures, contiguity and thresholds re-checked.
    The Sturm and containment code (m2_box_tools) is shared by the planner and the second check, so for that part
    only the matrix and the inputs are independent.
  Negative controls: a duplicated row (singular A), a box shifted off its branch, and a tampered entry of M all FAIL.
  Algebraic multiplicity ≥ 5 (the ray dispersion relation, item (M)): now also proved without SymPy.
    scripts/verify_m2_dispersion_fractions.py shows det(x − P) = (x² − s² − 1)⁴(x² − (1 + s²)/4)² Q/lc_x(Q) by exact
    evaluation on a 17 (s) × 65 (η) grid, every x-coefficient. After clearing det S⁰¹⁶ both sides have degree
    ≤ 16 in s and ≤ 64 in η, so this proves the identity. Q comes from scripts/m2_dispersion_export.py
    (data/m2_dispersion.json); a tampered Q fails.
  Still SymPy only: rotation equivariance (used only for the corollary about rotated covectors) and the general
    curved-symbol identity.
  A first second-check version with naive interval arithmetic could not confirm 46 boxes (too crude a bound; F-6).

Reproduce
  python3 -I scripts/m2_family_extended.py --out results/... --json data/m2_extended_certificate.json   (~20 s)
  python3 -I scripts/verify_m2_extended_fractions.py data/m2_extended_certificate.json --out results/...  (~10 s)
  Files, SHA-256 (run of 2026-11-22):
    b5ac2af05526dae644bf2bf3ca17c9b0d8e3e0e38cc7867d4593a2370dfe5cd9  scripts/m2_box_tools.py
    b8699e3417113d4021876a1a727e61bfb55ac14b9e244cf6a68f3e4e40b99999  scripts/m2_family_extended.py
    e95afcc71ed2baf84922d0dc86ec681ff914b3117092ff4f8a815e0e8c20c5c5  scripts/verify_m2_extended_fractions.py
    7bcf43c67673ddca843d4e8a4607d4753d1512885d46fbc6682c5bcf3a51f082  data/m2_extended_certificate.json
    49c32f869879b3d0c739522f07b755e9ec2334fc2dd78ccc655d230c0f92f3df  results/m2_family_extended_2026-11-22.txt
    17954800320e469e79dbcc4f16ab4ade5be1180067c1ad56c1d63e11563433a4  results/verify_m2_extended_fractions_2026-11-22.txt
    (dispersion: scripts/m2_dispersion_export.py, scripts/verify_m2_dispersion_fractions.py, data/m2_dispersion.json,
     results/m2_dispersion_export_2026-11-22.txt, results/verify_m2_dispersion_fractions_2026-11-22.txt; ~12 min on 2 cores)
  Python 3.13.18, SymPy 1.14.0, mpmath 1.3.0. Exit codes 0.

What this does NOT establish
  - The two small gaps near ε_L and ε_U; anything for |ε| ≥ ε_U, where det S⁰ also fails at 8.82.
  - Anything off this one-parameter curve (m, λ, b not varied); uniformity in ε; which property fails; physics.
  - Independent (non-AI) confirmation.

Review status: AI-checked only. Not peer reviewed, no human expert review, not independently replicated,
not formally verified.
```

## In plain words

Off the slice, the paper's crossing splits into a lower and an upper crossing. The lower one exists only while |ε| < 2.667, and its direction collapses onto the S²(3) axis at the end. The upper one exists while |ε| < 6.291, and its direction runs off to the S²(1) axis. The certificate shows there is a Jordan block at the lower crossing for all |ε| up to 2.66, and at the upper one for all |ε| up to 6.28. That covers almost all of each crossing's life. It is done with 108 small exact checks plus the earlier one near the slice.
