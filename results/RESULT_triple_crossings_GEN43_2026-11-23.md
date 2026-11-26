# Result: Jordan type J₃ ⊕ 3J₁ at the four real triple crossings, certified by the GEN43 engine   (date: 2026-11-23; recorded here 2026-11-25)

This is a pointwise result at four exact points, outside the split-crossing family. Unlike the other RESULT files, it was **certified by an external engine** (GEN43 r10.2.1), which is not part of this repository. Its certificate files are included in `data/gen43_r10_2_1_triple_crossings/`, but no script here re-derives them. See CONVENTIONS.md, "Results certified outside this repository".

```
Outcome type: POSITIVE (constrained Jordan structure certified), at each of four points separately
              Jordan type of the light eigenvalue on ker L: exactly J₃ ⊕ 3J₁
              Status in this repository: CERTIFIED EXTERNALLY

Parameter points (common: λ² = 1/2, m² = 823/142944, b₁ = 1, b₂ = 1/3, π′ = 1/2; only q changes)
  q = (1 + ε)·q_slice, q_slice = λ√(2/3)·m;  η = 35736·(q − q_slice), η² = 2450894·v, v = ε²
  v is a root of p₅, the squarefree part of Res_u(disc_X Q, C)/v² (degree 5, five real roots,
  isolated in (32,33), (54,55), (76,77), (102,103), (108,109)). The real triple crossings
  (u* > 0) are the roots
    v ≈ 102.677  (ε ≈ 10.133)   and   v ≈ 108.096  (ε ≈ 10.397),
  each with η > 0 and η < 0: four points. (The root near 32.29 is the exact exclusion at
  ε ≈ 5.68 in audit/RESULTANT_AND_ROOT_ISOLATION.md.)
  Exact or approximate: exact algebraic numbers (roots of p₅ with sign choices), isolated by
  Sturm sequences; ε values above are approximate labels only.

Crossing
  Covector: ray 1, ξ = (s, 0, 1, 0) with s² = u*(v) > 0 (real).
  Eigenvalue: the light eigenvalue r = |ξ| = √(1 + u*), real.
  Number field: K = ℚ(i)[t, η, s, r]/(p₅(t), η² − 2450894t, s² − u*(t), r² − 1 − u*(t)),
    degree 80 over ℚ, certified to be a field. Its 80 embeddings cover both triple
    crossings and both signs of η.

Constrained structure (on ker L(ξ))
  dim ker L(ξ) = 12
  Algebraic multiplicity of r on ker L: 6  (ray dispersion identity, quotient factor and
    model identities; same Q as data/m2_dispersion.json, sha256 baed27fd…)
  rank [L; P − r] = 12, so the eigenspace in ker L has dimension 4:
    ≥ 12 exact (prime witness in K, every embedding); ≤ 12 from the upper rank identity on ray 1
  rank [L; (P − r)²] ≥ 11, so dim ker (P − r)² on ker L ≤ 5:
    exact (prime witness) and, separately, rigorous intervals
  With alg 6 and geo 4 the type is J₃ ⊕ 3J₁ (kernel of (P − r)² has dimension 5) or
    2J₂ ⊕ 2J₁ (dimension 6). The rank bound excludes the second. Result: J₃ ⊕ 3J₁.
  Chain: not written out explicitly. The certificates prove the ranks, which imply a
    chain v → w → u inside ker L exists.   files: data/gen43_r10_2_1_triple_crossings/
  Engines: GEN43 r10.2.1 (its own exact Q(i)/number-field arithmetic and rigorous intervals).
  Cross-checked in a second engine: within GEN43, two routes (exact and interval) for
    the classification. NOT re-derived by any script in this repository.

Hypotheses re-checked at this point (the paper's argument needs all four)
  det S⁰ ≠ 0: CERTIFIED EXTERNALLY at each point (GEN43, rigorous interval enclosure of
    width about 1e-75). This matters here: det S⁰ vanishes at ε ≈ 10.118, close to the
    triple crossing at ε ≈ 10.133. (step2 log: det S⁰ as an exact function of ε.)
  Reduction identities with q free: covered by the general statement in this repository
    (results/step2_reduction_q_free_2026-11-22.txt: symbolic q, valid wherever det S⁰ ≠ 0).
    Not re-run at these points.
  Coercivity on some ball: covered by the general statement in this repository
    (audit/OFFSLICE_COERCIVITY_FINAL.md, Proposition C and Theorem T: witness m, λ, b,
    every real q with det S⁰(q) ≠ 0). Not re-checked at these points.
  Curved-symbol identification: covered by the general statement in this repository
    (results/m2_curved_symbol_general_2026-11-22.txt: every covector and every q, SymPy
    only). Not re-checked at these points.

Conditional consequence (not claimed as a result)
  Theorem T applies to any real q with det S⁰(q) ≠ 0 and a Jordan chain of P on ker L at a
  real eigenvalue. If the external certificate is accepted, its hypotheses hold at these four
  points, so the loss-free package WP^ℓ(R, c, T) of the parent paper would fail there too.
  This file records that implication only conditionally, because one input is certified
  outside this repository (CONVENTIONS.md).

Reproduce
  Not reproducible from this repository alone (the engine is not published here).
  Engine: GEN43_engine_r10_2_1_2026-11-23_1957.zip,
    sha256 56f0a3a8cce07dd5848672f287947f486f09b94ba74ec2c6a95ec2d52e36b5b2
  Run: by the owner on his own PC (Ryzen 9 9900X), 2026-11-23, profile r10
    (157 results, 0 mathematical failures, 0 software errors).
  Anchor: SHA-256 of the full results-folder manifest
    717c130f94bf8174c39886d3ddfac0a99c6340d2c4d10db59457de63b9c1524a
    (stored by the owner outside the folder; recomputed 2026-11-25: MATCH).
  Files: the 34 files in data/gen43_r10_2_1_triple_crossings/. They are byte-identical
    to the run (each checked against the run's certificate_index.json on 2026-11-25;
    see CERTIFICATE_INDEX_SUBSET.json).
  Stricter re-check: GEN43 r10.3 Gate 1b (2026-11-25), strict gate v2. All four
    classifications remain VERIFIED_EXACT, with 0 label changes. The scope entries they use
    were approved by the owner. That approval is a scope judgement on the records, not a
    review of the mathematics.
  Shared inputs: GEN43 transcribes the defining equations from the parent paper and code,
    and uses a different gamma representation and its own arithmetic. Agreement with this
    repository is therefore "shared formulas, different implementation", not independent
    confirmation.

What this does NOT establish
  - Anything at other points, or on an open set of parameters.
  - Anything for π′ ≠ 1/2 (the algebraic-multiplicity input uses π′ = 1/2).
  - An explicit chain written out in this repository.
  - Uniformity, or which of existence, uniqueness or the bound fails.
  - Anything about physics.
  - Independent or human verification.

Review status: AI-checked only. Not peer reviewed, no human expert review, not formally
verified, not independently replicated. Supersedes, as far as the Jordan type is concerned,
the EXPLORATORY triple-crossing entry in audit/README.md (that entry is kept unchanged).
```

A size-3 Jordan block is a stronger algebraic defect than the J₂ found on the split-crossing family. Whether it changes anything beyond the conditional consequence above, for example the severity question in the parent ROADMAP §3, has not been worked out.

GEN43 r10.2.1 also certified the light-factor chain as VERIFIED_EXACT, for every parameter point with det S⁰ ≠ 0, every π′ and every ξ with ξ·ξ ≠ 0:
- L·Y = 0;
- L·P₀ = C₁·L;
- the constrained light factor.

That is a wider statement than this repository's grid proof at π′ = 1/2 (`audit/RESTRICTED_SYMBOL_COMPLETE_PROOF.md`). Only a summary is recorded here; no files are included.
