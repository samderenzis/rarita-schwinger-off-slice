# RESTRICTED_SYMBOL_COMPLETE_PROOF — invariance of ker L, characteristic polynomial, multiplicities

Date: 2026-11-22 (GEN43 Phase 1 final form; supersedes `RESTRICTED_SYMBOL_PROOF.md`, kept unchanged). Author: Claude (AI). Review status: AI-checked only. Tags as in `OFFSLICE_COERCIVITY_FINAL.md`.

## Question
Make explicit every step of: invariant constraint space → restricted characteristic polynomial → algebraic multiplicity → geometric multiplicity → Jordan block. Account for elimination factors, denominators, spurious roots and the quotient part.

## Conventions
- Witness m, λ, b₁ = 1, b₂ = 1/3; q = q_slice + η/35736 with η = ε√2450894 (q = (1 + ε)q_slice). d₀(η) := det S⁰ ≠ 0 is assumed throughout (true except at ε ≈ ±8.823 and ±10.118).
- P(ξ) = M⁻¹B(ξ) with the paper's M = [Γ^{0ik}] (code `MMmat`; code `Mmat` = −M). P is 16×16, linear in ξ.
- L(ξ) = [Σ_j ξ_jΓ^{0jk}]_k is 4×16, linear in ξ, and q-free. C₁(ξ) = π′Σ_a ξ_aΓ⁰Γ^a is 4×4.
- Ray: ξ = (s, 0, 1, 0), |ξ|² = 1 + s². Parametrisation s = 2τ/(1 − τ²), |ξ| = (1 + τ²)/(1 − τ²). u = s², v = ε².

## Summary theorem (R)
For the witness m, λ, b, any π′, any ξ ≠ 0 on the ray and any ε with det S⁰ ≠ 0:
(i) ker L(ξ) is 12-dimensional and P(ξ)-invariant; (ii) det(x − P) = det(x − A)(x² − π′²|ξ|²)²; (iii) det(x − A) = (x² − |ξ|²)⁴ Q/lc_x(Q), independent of π′; (iv) alg_A(|ξ|) = 4 + mult_X(|ξ|²; Q); (v) geo_A(|ξ|) ≥ 4 everywhere, with equality on the certified boxes. Together with the resultant analysis (`RESULTANT_AND_ROOT_ISOLATION.md`) these give the Jordan classification (`JORDAN_STRUCTURE_COMPLETE.md`).
**Status:** (i)–(iv) ESTABLISHED and REPRODUCED (π′ = 1/2 in the second engine; other π′ via R5, SymPy only). (v) ESTABLISHED; the upper bound on rank is proved in the fractions engine; the box bounds have a second checker without SymPy that shares the containment routine.

## Steps

**R1. rank L(ξ) = 4 and dim ker L(ξ) = 12 for every ξ ≠ 0.** From L L† = 3|ξ|²I₄, L is surjective.
ESTABLISHED (symbolic ξ, SymPy) and REPRODUCED (fractions engine, along the ray and at four points).

**R2. Invariance: L(ξ)P(ξ) = C₁(ξ)L(ξ), hence P(ξ) ker L(ξ) ⊆ ker L(ξ).** If Lv = 0 then LPv = C₁Lv = 0.
- SymPy (step 2): symbolic ξ, π′ and q at the witness m, λ, b.
- Fractions engine (`audit/restricted_symbol_fractions.py` (I)): at π′ = 1/2, for all ξ ∈ ℝ⁴ and all η. d₀(LP − C₁L) has degree ≤ 2 in each ξ_j and ≤ 4 in η, and is checked on a 3⁴ × 5 grid with d₀ ≠ 0.

ESTABLISHED, REPRODUCED (π′ = 1/2).

**R3. Block form and the quotient.** Let V be any complement of ker L. L|_V : V → ℂ⁴ is an isomorphism (R1). In a basis adapted to ker L ⊕ V,
  P = [[A, B], [0, D]],  A = P|ker L (12×12),
and L(Pv) = C₁(Lv) shows D ≅ C₁ under L|_V. Hence
  det(x − P) = det(x − A)·det(x − C₁).
ESTABLISHED (linear algebra).

**R4. det(x − C₁) = (x² − π′²|ξ|²)².** C₁² = π′²|ξ|²I and tr C₁ = 0 (paper Prop. 5.1(c); fractions engine (II)). So the eigenvalues are ±π′|ξ|, each twice; if π′ = 0, C₁ = 0.
ESTABLISHED, REPRODUCED.

**R5. Independence of π′.** P_{π′} − P_{π″} = (π′ − π″)M⁻¹TL vanishes on ker L, so A does not depend on π′.
ESTABLISHED (step 2, SymPy; paper Prop. 5.1(a)).

**R6. Restricted characteristic polynomial along the ray.** At π′ = 1/2, for all s and all η with d₀ ≠ 0,
  det(x − P) = (x² − s² − 1)⁴ (x² − (1 + s²)/4)² Q(x, s, ε)/lc_x(Q),
with lc_x(Q) = (9803576ε² − 763140625)(39214304ε² − 4014616321), a constant multiple of d₀.
- SymPy: ker L route, with ε symbolic (`scripts/m2_dispersion_export.py`).
- Fractions engine: the full 16×16 characteristic polynomial on a 17 × 65 grid with degree bounds ≤ 16 in s and ≤ 64 in η after clearing d₀¹⁶ (`scripts/verify_m2_dispersion_fractions.py`; derivation in DEGREE_BOUND_DERIVATION, to be written).

Dividing by det(x − C₁) = (x² − (1 + s²)/4)² (R3, R4 with |ξ|² = 1 + s²):
  **det(x − A) = (x² − s² − 1)⁴ Q(x, s, ε)/lc_x(Q).**
Both sides are monic of degree 12 in x. No elimination or denominator factor remains: the only denominator, lc_x(Q), is the nonzero d₀ up to a constant. By R5 this holds for every π′.
ESTABLISHED, REPRODUCED.

**R7. Algebraic multiplicity at r = |ξ| > 0.** Q is a polynomial in X = x², Q = aX² + bX + c (checked). So mult_x(r; Q) = mult_X(r²; Q), and
  alg_A(r) = 4 + mult_X(1 + s²; Q).
Also r² − r²/4 ≠ 0, so the quotient factor does not vanish at r, and alg_P(r) = alg_A(r).
- A crossing means C(u, v) := Q(X = 1 + u) = 0. Then alg_A(r) ≥ 5.
- alg_A(r) = 5 exactly iff disc_X Q(u, v) ≠ 0. Otherwise both sheets meet the light cone and alg_A(r) = 6.

ESTABLISHED.

**R8. Geometric multiplicity.** K(s) (entries a + bs) satisfies L(s,0,1,0)K(s) = 0 for all s and contains a 12×12 identity submatrix (rows 4–15). So it is a basis of ker L for every s (fractions engine (III)). Hence
  geo_A(r) = 12 − rank((P − r)K) = 12 − rank M(τ, η),  M := D_n(P − |ξ|)K,
where D_n = (τ² − 1)²(4η² − 27625²)(16η² − 63361²) ≠ 0 for τ ∈ (−1, 1) and |η| < 13812.5.
- **Upper bound rank M ≤ 8 everywhere:** every 9×9 minor of M has degree ≤ 36 in τ and in η, because each entry has degree ≤ 4 in each variable; this entry bound is proved in `scripts/verify_m2_extended_fractions.py`. The minors vanish on a 37 × 37 grid (fractions engine (IV)), hence identically. So geo_A(|ξ|) ≥ 4 at every point of the ray.
- **Lower bound rank M ≥ 8:** on the certified boxes. Ray 1 (ξ = (s,0,1,0)): `RESULT_m2_q_direction`, `RESULT_m2_extended` (108 boxes, |ε| ≤ 5.7 upper / 2.66 lower) and the gap boxes of `EPSILON_ENDPOINT_CLOSURE.md` (lower branch up to ε_L, end box containing s = 0). Ray 2 (ξ′ = (1,0,t,0), used for the upper branch near ε_U, end box containing t = 0): hand-pivoted basis K′ of ker L(ξ′), valid for every real t including 0 because the pivot block has det = (1 + t²)². Each box bound is ‖CA − I‖_∞ < 1 for an 8×8 submatrix A of M with a rational approximate inverse C; on the boxes the largest bound is 0.9708 (extended) and 0.1635 (gaps).

ESTABLISHED (upper bound REPRODUCED in the sense of the fractions engine only; box bounds have a second check without SymPy).

## Consequences
1. **Off the crossings the light eigenvalue is semisimple.** At any (s, ε) on the ray with d₀ ≠ 0 and Q(1 + s²) ≠ 0: alg_A(|ξ|) = 4 (R7) and geo_A(|ξ|) ≥ 4 (R8), so geo = alg = 4. ESTABLISHED.
2. **At a crossing** with disc_X Q ≠ 0 and rank M = 8: alg = 5 and geo = 4, so there is exactly one Jordan block of size 2 (see JORDAN_STRUCTURE_COMPLETE.md).
3. By rotation equivariance (SymPy only) and P(−ξ) = −P(ξ), the same holds on the whole crossing cones {|ξ_a|²/|ξ_b|² = u±(ε)} and at the eigenvalue −|ξ|.

## Remaining gaps
- R5 and the rotation corollary are SymPy only.
- R2 in the fractions engine covers π′ = 1/2 only. Other π′ rely on R5.
- The degree-bound derivations are summarised here and in the scripts' docstrings; a separate DEGREE_BOUND_DERIVATION.md is planned (Phase 6).
- Every crossing covector lies on the SO(2)×SO(2) orbit of the ray (the crossing cones are the sets |ξ_a|²/|ξ_b|² = u±(ε)), so the ray results cover all crossings once rotation equivariance is granted; that equivariance is SymPy only.
