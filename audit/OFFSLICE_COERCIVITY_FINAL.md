# OFFSLICE_COERCIVITY_FINAL — the analytic layer off the slice: coercivity, necessity lemma, reduction, main theorem

Date: 2026-11-22 (GEN43 Phase 1 final form; supersedes `OFFSLICE_COERCIVITY_AUDIT.md`, which is kept unchanged as the earlier draft). Author: Claude (AI). Review status: AI-checked only. Not peer reviewed, no human expert review, not independently replicated.

## Status tags used
| Tag | Meaning in this project |
|---|---|
| ESTABLISHED | proved by a written argument and/or an exact computation in this project |
| REPRODUCED | ESTABLISHED, and also recomputed exactly by a second, separately written engine of this project. This is still the same researcher and AI process, not independent replication. |
| SUPPORTED | argued in writing, but some step relies on an unchecked reading or a citation |
| EXPLORATORY | numerical or heuristic only |
| INCOMPLETE | a known gap remains |
| DISPROVED | shown false |

Paper references are to paper version 2 (main repository v1.0.11, `paper/main.tex`). The paper has no section numbers in its statements, so the propositions are named here: "Prop. 4.1" = Proposition *Estimate on constrained test functions* (`prop:wh`), "Prop. 4.2" = Proposition *Small-ball coercivity* (`prop:coer`).

## 1. Question
Fix the background (b₁ = 1, b₂ = 1/3, chart V′), m > 0, λ ∈ ℝ, π′, and an arbitrary real charge q (not necessarily q = λ√(2/3)m). Which steps of the paper's proof chain still hold, with which constants, and what extra hypotheses are needed?

## 2. Where the slice relation is used
A complete reading of `main.tex` (§2–§5, Appendices A–B) finds **one** use of q = λ√(2/3)m: the parameter choice (eq. `witness-params`), which makes the witness data lie in ℚ(i). No lemma, proposition or proof step uses it. The table lists every step and what it needs from q.

| Step (paper) | Uses the slice relation? | What it needs, for a given q | q-dependent quantities | Status off the slice |
|---|---|---|---|---|
| Field equations, covariant derivative (§2.2) | no | nothing (q is a free real parameter in the definitions) | iqA_B in ∇ | ESTABLISHED (definition) |
| Lemma INT (§3.1) | no | nothing: X_A, Ŷ_A contain λ, m, F only | none | ESTABLISHED; X, Ŷ checked q-free (step 2) |
| Divergence identity, Prop. 3.2 | no | the written proof is for any q: the term (iq/2)Γ^{ABρ}(dA)_{AB} comes from the curvature of the U(1) connection | S^ρ affine in q | ESTABLISHED (written proof); exact check with q symbolic on the background (step 2, SymPy) |
| Elimination of ψ₀, Lemma 3.3 | no | **det S⁰(q) ≠ 0** | S⁰, 𝒲 = −(S⁰)⁻¹S^k | ESTABLISHED for every q with det S⁰ ≠ 0. Witness m, λ: det S⁰ = (9803576ε² − 763140625)(39214304ε² − 4014616321)/132101812606622109696, which vanishes only at ε² = 27625²/9803576 (ε ≈ 8.823) and ε² = 4014616321/39214304 (ε ≈ 10.118) |
| Constraint intertwining, Prop. 3.4 | no | written proof for any q; exact check with q symbolic (step 2) | 𝒞, ⋆ contain iqA_a | ESTABLISHED |
| Equivalence and transfer lemmas (§3.5) | no | det S⁰ ≠ 0; 𝒲 constant (S^ρ are constant frame matrices for every q, because ∇X = 0 and dA = F) | none beyond 𝒲 | ESTABLISHED |
| Prop. 4.1 (estimate on constrained test functions) | no | the package with constants K, α; Prop. 3.4; smooth coefficients of K′ | none: the constants come from the hypothetical package | ESTABLISHED |
| Prop. 4.2 (small-ball coercivity) | no | the bound m₀(q) = sup‖M₀‖ | **only m₀ and hence r₀** (§3 below) | ESTABLISHED with explicit constants (§3, §4) |
| Necessity lemma hypotheses (Pos), (BT) | no | μ and the cone geometry | none | ESTABLISHED |
| (Reg_s) | no | μ, L_j, L₀ smooth on V̄′; A_j Lipschitz and 𝔅 bounded on B̄ | L₀ affine in q; A_j, 𝔅 contain (S⁰)⁻¹S^k and iqA | ESTABLISHED for every q with det S⁰ ≠ 0 (all coefficients are smooth on V̄′, where sin θᵢ ≥ 1/2) |
| (C1) | no | Prop. 4.2 | r₀(q) | ESTABLISHED (§4) |
| (H_X) | no | Prop. 4.1; p_Λ ∈ H^{ℓ+2} (Lemma *Higher interior regularity* with s = ℓ + 2, which needs (Reg_s) for that s) | none | ESTABLISHED |
| Lemma *Oscillatory H⁻¹ bound* | no | weight μ only | none | ESTABLISHED |
| Appendix B (Lax–Milgram, adjoint identity, H², Hˢ, correction bounds, proof) | no | (Pos), (Reg_s), (C1), (BT), (H_X); constants N₁, N₂ depend on M_j, M₀ | through M₀ (finite for each q) | ESTABLISHED for each fixed q |
| Curved-symbol identification at y₀ | no | the principal symbol of the curved K′ at y₀ equals iP(ξ), ξ = E(y₀)ᵀk | P | ESTABLISHED for every k and every q at the witness m, λ, b (`scripts/m2_curved_symbol_general.py`, SymPy only; not REPRODUCED) |
| Proof of the Main Theorem (§5.3) | no | a chain v → w in ker L(ξ) at a real eigenvalue r, at the same y₀ and k | v, w, r | ESTABLISHED given the chain |

Nothing in the chain requires a spectral gap, a symmetriser, strong hyperbolicity or a frequency estimate beyond the oscillatory H⁻¹ bound. The invariance P ker L ⊂ ker L is not used by the necessity lemma (the paper states this explicitly).

## 3. Exact q-dependence of the coercivity constants
Computed by `audit/coercivity_constants.py` (SymPy, exact; m, λ, q symbolic; log `results/audit_coercivity_constants_2026-11-22.txt`).

- R⁰Φ = Σ_j L_j(y)∂_jΦ + L₀(y)Φ exactly. R⁰ contains neither ψ₀ nor ∂_t. **ESTABLISHED**
- L_j(y) = Σ_a E_a^j(y)[Γ^{0ak}]_k contains **no q, m or λ**. **ESTABLISHED**
- L(ξ)L(ξ)† = 3|ξ|²I₄ for symbolic ξ; it contains no q. **ESTABLISHED, REPRODUCED** (also at four points and along the ray in the fractions engine)
- M₀ := L₀† − Σ_j μ⁻¹∂_j(μL_j†) depends on y only through cᵢ = cot θᵢ. It is **affine in q**, M₀ = M₀⁽⁰⁾ + qM₀⁽¹⁾, and q enters only through the gauge term iqA_j (A = −cot θ₁e² − cot θ₂e⁴). **ESTABLISHED**
- Exact Frobenius identity (the terms linear in q cancel):

  ‖M₀(y; q)‖²_F = 3c₁² + c₂²/3 + 10λ²/3 + 16m² + 12q²(c₁² + c₂²).   **ESTABLISHED (SymPy only)**

- On V̄′, θᵢ ∈ [π/6, 5π/6], so |cᵢ| ≤ cot(π/6) = √3. Hence, for **every** real q, m, λ:

  m₀(q) := sup_{V̄′}‖M₀‖₂ ≤ sup_{V̄′}‖M₀‖_F = **√(10 + 10λ²/3 + 16m² + 72q²)**.

  At the witness m, λ: m₀(q) ≤ √(5746777632q² + 938543502)/8934. On the slice this is ≈ 3.449.

- q-free constants: σ = √3·min(b₁, b₂) = 1/√3; ℓ₁ = 6√(b₁² + b₂²) = 2√10; μ₋ = inf μ = 1/(4b₁²b₂²) = 9/4; Poincaré constant 1/j₁,₁ < 1.

## 4. Standalone off-slice coercivity proposition

**Proposition C (coercivity for arbitrary real charge).** Let b₁ = 1, b₂ = 1/3, m > 0, λ ∈ ℝ, q ∈ ℝ, and put

  m̄(q) := √(10 + 10λ²/3 + 16m² + 72q²),  r₀(q) := min(π/3, σ/(2(ℓ₁ + m̄(q)))),  σ = 1/√3, ℓ₁ = 2√10,
  c₀ := (σ/2)·μ₋^{1/2}·(1 + π²/9)^{−1/2} = (√3/4)(1 + π²/9)^{−1/2} ≈ 0.2990.

Then for every ρ′ ≤ π/6 with 2ρ′ ≤ r₀(q), the geodesic ball B′ = B_d(y₀, ρ′) satisfies

  ‖R^{0*}ρ‖_{L²(B′, μ)} ≥ c₀‖ρ‖_{H¹}  for all ρ ∈ H¹₀(B′).

c₀ does not depend on q. r₀ depends on q only through m̄(q). It is uniform on compact q-sets and bounded below by a constant times 1/|q| as |q| → ∞.

*Proof.* (Following the paper's steps (a)–(f), with every constant made explicit.)
1. **Symbol bound.** For ξ = E(y)ᵀk, |ξ| ≥ min(b₁, b₂)|k| because 1/sin θᵢ ≥ 1. From L(ξ)L(ξ)† = 3|ξ|²I we get |L(y,k)†v| = √3|ξ||v| ≥ σ|k||v|.
2. **Frozen coefficients.** On the Euclidean chart ball 𝒰_r = B_e(y₀, r), r ≤ π/3, Plancherel gives ‖Σ_j L_j(y₀)†∂_jρ‖_{L²(dy)} ≥ σ‖∇ρ‖ for ρ ∈ C^∞_c(𝒰_r), and by density for ρ ∈ H¹₀(𝒰_r).
3. **Freezing error.** Only E₂^{φ₁} = b₁/sin θ₁ and E₄^{φ₂} = b₂/sin θ₂ vary. Since |∂_θ(1/sin θ)| = |cos θ|/sin²θ ≤ 2√3 on [π/6, 5π/6], we have |E(y) − E(y₀)| ≤ 2√3 b·r on 𝒰_r. Each [Γ^{0ak}]_k has operator norm ≤ √3 (three unitary blocks), so (Σ_j‖L_j(y) − L_j(y₀)‖²)^{1/2} ≤ 6r√(b₁² + b₂²) = ℓ₁r. By Cauchy–Schwarz, ‖Σ_j(L_j(y) − L_j(y₀))†∂_jρ‖ ≤ ℓ₁r‖∇ρ‖.
4. **Lower order.** ‖M₀ρ‖ ≤ m₀(q)‖ρ‖ ≤ m̄(q)‖ρ‖ (§3). Poincaré on a ball of radius r in ℝ⁴ (first Dirichlet eigenvalue (j₁,₁/r)², j₁,₁ ≈ 3.83 > 1) gives ‖ρ‖ ≤ r‖∇ρ‖.
5. **Absorption.** R^{0*}ρ = −Σ_j L_j†∂_jρ + M₀ρ, so ‖R^{0*}ρ‖_{L²(dy)} ≥ (σ − (ℓ₁ + m̄(q))r)‖∇ρ‖ ≥ (σ/2)‖∇ρ‖ for r ≤ r₀(q). Then ‖·‖_μ ≥ μ₋^{1/2}‖·‖_{L²(dy)} and ‖ρ‖²_{H¹} ≤ (1 + r²)‖∇ρ‖² ≤ (1 + π²/9)‖∇ρ‖² give the claim with c₀.
6. **Geodesic balls.** The chart metric diag(1/b₁², sin²θ₁/b₁², 1/b₂², sin²θ₂/b₂²) is ≥ ¼I on V′, so B_d(y₀, ρ′) ⊂ B_e(y₀, 2ρ′) ⊂ V′ for ρ′ ≤ π/6. Extension by zero maps H¹₀ of the geodesic ball into H¹₀(𝒰_{2ρ′}). ∎

**Status: ESTABLISHED** (written proof; the constant m̄(q) from an exact SymPy computation, not REPRODUCED in the fractions engine).

Explicit values at the witness m, λ (lower bounds for r₀): 0.0295 on the slice; 0.0294 at ε = 1; 0.0270 at ε = 6.28; 0.0187 at q = 1; 0.0032 at q = 10.

## 5. The other analytic hypotheses for arbitrary q
**Proposition N (necessity lemma hypotheses off the slice).** Fix q with det S⁰(q) ≠ 0, and the radii r*(q) := min(π/6, r₀(q)/2), B′ = B_d(y₀, r*), B″ = B_d(y₀, 3r*/4), B = B_d(y₀, R) with R < r*/2, and cones B_t as in the paper. Then (Pos), (Reg_s) for every s, (C1), (BT) hold. If WP^ℓ(R,c,T) holds, then (H_X) holds on X* with K₀ = K.

*Proof.* (Pos): μ = sin θ₁ sin θ₂/(b₁²b₂²) ∈ [9/4, 9] on V̄′. (Reg_s): μ, L_j are smooth and q-free; L₀ is a polynomial in cot θᵢ with coefficients affine in q; A_j and 𝔅 are built from E(y), the constant matrices (S⁰)⁻¹S^k, M⁻¹, and the smooth coefficients of R; all are smooth on V̄′ when det S⁰ ≠ 0. (C1): Proposition C. (BT): the cone geometry does not involve q. (H_X): Prop. 4.1, whose proof uses only the package, Prop. 3.4 (valid for every q) and the smoothness of the coefficients of K′. The wave packets lie in X* by the Hˢ interior regularity lemma with s = ℓ + 2. ∎   **Status: ESTABLISHED**

## 6. Off-slice main theorem (transfer)

**Theorem T.** Let b₁ = 1, b₂ = 1/3, λ² = 1/2, m² = 823/142944, any π′, and a real q with det S⁰(q) ≠ 0. Suppose that for some ξ ≠ 0 and real r there are v, w ∈ ℂ¹⁶ with Lv = Lw = 0, (P(ξ) − r)v = 0, (P(ξ) − r)w = v, v ≠ 0. Then for every R ∈ (0, r*(q)/2), every integer ℓ ≥ 2, every c > 0 and every T ∈ (0, R/c), the package WP^ℓ(R,c,T) is false. Here r*(q) = min(π/6, r₀(q)/2) is given explicitly by Proposition C.

*Proof.* The paper's proof of the Main Theorem, with the witness replaced by (ξ, r, v, w), and with the necessity lemma applied at y₁ = y₀ and k = E(y₀)^{−T}ξ. Its hypotheses hold by Proposition N. Its symbol at (y₀, k) is P(ξ) by the curved-symbol identification for every k and q (`scripts/m2_curved_symbol_general.py`). ∎

**Status: ESTABLISHED** for the witness m, λ, b and every q with det S⁰ ≠ 0. The input with the weakest verification is the general curved-symbol identity (SymPy only). For other m or λ, the curved-symbol identity has not been checked (INCOMPLETE).

## 7. Appendix B constants off the slice
The paper's Appendix B constants N₁, N₂ (adjoint identity, H² and Hˢ bounds, correction bounds) enter only through sup‖M₀‖ and sup‖∂M₀‖ on V̄′, besides q-free geometric quantities. `audit/coercivity_constants.py` section (6) (SymPy, exact) shows:
- M₀ is affine in c₁ and c₂ separately (no c₁c₂ and no cᵢ² terms). **ESTABLISHED (SymPy only)**
- ‖∂M₀/∂c₁‖²_F = 3 + 12q² and ‖∂M₀/∂c₂‖²_F = 1/3 + 12q². **ESTABLISHED (SymPy only)**
- ∂/∂θᵢ = −csc²θᵢ ∂/∂cᵢ with csc² ≤ 4 on V̄′, and ∂_φ M₀ = 0. Hence ‖∂_{θ₁}M₀‖ ≤ 4√(3 + 12q²), ‖∂_{θ₂}M₀‖ ≤ 4√(1/3 + 12q²).

So N₁ ≤ (q-free) + m̄(q) and N₂ ≤ (q-free) + 4√(3 + 12q²) + 4√(1/3 + 12q²). Both are finite for every q, uniform on compact q-sets, and grow at most linearly in |q|. **ESTABLISHED** (bounds), with the q-free parts as in the paper.

## 8. Pointwise, compact-set and full-axis statements
| Constant | Pointwise in q | Uniform on compact q-sets | Uniform on all of ℝ |
|---|---|---|---|
| c₀ (coercivity) | yes | yes | **yes** (q-independent, √3/4·(1 + π²/9)^{−1/2}) |
| m̄(q), r₀(q), r*(q) | explicit, finite / positive | yes | **no** as bounds: r₀ ≳ 1/|q| as |q| → ∞ |
| N₁, N₂ (Appendix B) | finite | yes | **no** as bounds (linear growth) |
| (Reg_s), 𝒲 = −(S⁰)⁻¹S^k | needs det S⁰(q) ≠ 0 | on compact sets avoiding ε ≈ ±8.823, ±10.118 | **no**: blows up at the zeros of det S⁰ |
| K₀ (necessity lemma) | = package constant K | — | — (no q-dependence claimed) |

Every theorem in this project is a statement for each fixed q (pointwise). Uniform versions on compact q-sets follow from the table, but none is used. Nothing uniform on all of ℝ is claimed except c₀.

## 9. Remaining gaps
- The general curved-symbol identity is SymPy only (INCOMPLETE as a REPRODUCED item).
- The M₀ formula is SymPy only (ESTABLISHED, not REPRODUCED).
- For m, λ, b other than the witness values, Theorem T needs the curved-symbol identity at those values (not done).
- The written proofs here are AI-authored and AI-checked only.

## 10. Strongest justified conclusion
For the witness m, λ, b and every real q with det S⁰(q) ≠ 0, the paper's analytic framework applies unchanged, with an explicit, q-dependent radius and a q-independent coercivity constant. Any real Jordan chain inside ker L at such a q therefore rules out the loss-free constrained L² package on small cones around y₀. No uniformity in q is claimed. No physical, causality or hyperbolicity claim follows from this.
