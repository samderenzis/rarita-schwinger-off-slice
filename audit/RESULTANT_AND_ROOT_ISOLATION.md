# RESULTANT_AND_ROOT_ISOLATION — where the light crossing can meet a G-sheet coincidence

Date: 2026-11-22 (GEN43 Phase 1). Author: Claude (AI). Review status: AI-checked only; no human review; not independently replicated. Tags as in `OFFSLICE_COERCIVITY_FINAL.md`.

Script: `audit/jordan_classification.py` (SymPy 1.14, exact rational arithmetic). Log: `results/audit_jordan_classification_2026-11-22.txt` (script sha256 6cea97cd…, input `data/m2_dispersion.json` sha256 baed27fd…).

## 1. Question
At a crossing (Q(X = 1 + u) = 0, u = s² > 0), the light eigenvalue r = |ξ| has algebraic multiplicity 5 exactly when disc_X Q ≠ 0 (RESTRICTED_SYMBOL_COMPLETE_PROOF R7). For which ε² = v does a real crossing coincide with a zero of disc_X Q?

## 2. The two polynomials
Q is a quadratic in X = x² with coefficients in ℤ[u, v] (checked: Q is a polynomial in X, u, v).
- disc_X Q = 3129942848356224 · v · D(u, v), with D of degree 2 in u and 3 in v (explicit in the log). The factor v is the slice degeneracy (both G-sheets coincide at ε = 0).
- C(u, v) := Q(X = 1 + u) = 8868484 · (54863649u²v² − 7763306328u²v + 221321320704u² + 12191922uv² − 1372250448uv − 18065203200u + 677329v² − 56628984v + 368640000).

## 3. Resultant and factorisation
**Lemma 1.** Res_u(disc_X Q, C) = κ · v² · p₅(v)² for a nonzero rational κ, where

p₅(v) = 377571474600343v⁵ − 141153925904173916v⁴ + 20332537164974086802v³ − 1400990293471261118908v² + 45801950198017653991053v − 562553552914350795281250.

*Proof.* Exact SymPy computation of the resultant and comparison with κv²p₅² (identity checked coefficient by coefficient). ∎ **ESTABLISHED (SymPy only).**

Since the resultant lies in the ideal (disc_X Q, C) ⊂ ℚ[u, v], every common zero (u, v) ∈ ℂ² satisfies v = 0 or p₅(v) = 0. This needs no condition on leading coefficients: Res = a·disc + b·C with a, b ∈ ℚ[u, v].

**Lemma 2 (irreducibility).** p₅ is irreducible over ℚ.
*Certificate.* p₅ is primitive (content 1). Its leading coefficient is 823⁵ = 377571474600343, not divisible by 11. The reduction p₅ mod 11 is irreducible of degree 5 over 𝔽₁₁ (SymPy factorisation over GF(11)). A factorisation over ℚ would, by Gauss's lemma, give one over ℤ with factor degrees preserved mod 11. ∎ **ESTABLISHED** (the mod-11 check is a small finite computation that can be repeated by hand or in any CAS).

## 4. Real root isolation of p₅
Exact Sturm/Vincent isolation gives five real roots, one in each of
[32, 33], [54, 55], [76, 77], [102, 103], [108, 109],
that is ε ≈ 5.6826, 7.3816, 8.7347, 10.1330, 10.3969. **ESTABLISHED.**

Relative to the existence thresholds of the crossings (exact, from the u-discriminant of C):
- ε_L² = (34404 − 12√5659689)/823 ≈ 7.115 < 32: **no root of p₅ in (0, ε_L²].**
- ε_U² = (58228 − 20√1645345)/823 ≈ 39.58, and 32 < ε_U² < 54: **exactly one root v₀ in (0, ε_U²].**

Refined: v₀ ∈ I₀ := [187607633/5809776, 11905019/368671] (ε₀ ≈ 5.68257999307).

## 5. Branch identification at v₀
On I₀, C(·, v) has two real roots in u:
- The upper crossing u₊(v) ∈ [17/10, 19/10] for every v ∈ I₀. Proof by exact Sturm counts: exactly one root in [17/10, 19/10] at the left end of I₀; no root of C on the horizontal sides u = 17/10, 19/10 for v ∈ I₀ (Sturm in v); no root above 19/10.
- The other real root is the continued lower branch, in [−2/100, −1/100]. It is negative, so it is not a real covector.

On the box [17/10, 19/10] × I₀, an exact interval enclosure of disc_X Q is [3.419·10³⁶, 9.527·10³⁶], which excludes 0. So the common zero of disc and C at v₀ is the lower branch at u ≈ −0.0150. **ESTABLISHED.**

## 6. Conclusion
**Proposition RI.** For every real crossing (u > 0) with 0 < v ≤ ε_L² on the lower branch, or 0 < v ≤ ε_U² on the upper branch (endpoints included), disc_X Q(u, v) ≠ 0. Hence the light eigenvalue has algebraic multiplicity exactly 5 on ker L there.

At the endpoints themselves the crossing sits at s = 0 (lower) or s = ∞ (upper, i.e. ξ = (1,0,0,0) in ray-2 coordinates). Multiplicity 5 there is shown separately by `audit/endpoint_multiplicity.py` (homogenised Q; X = 1 is a simple root at both endpoints, with derivatives −67603406660674704 − 27990756539088√5659689 and 422979956239154544 + 340756598051280√1645345). See `EPSILON_ENDPOINT_CLOSURE.md`. **ESTABLISHED (SymPy only).**

## 7. Exploratory context (outside the existence ranges)
Common zeros of disc_X Q and C at the other roots of p₅ (high-precision numerics of the exact roots; status EXPLORATORY for anything not used above):
| ε | crossing roots u | common with disc | real covector? |
|---|---|---|---|
| 5.6826 | −0.0150, 1.7990 | −0.0150 | no |
| 7.3816 | −1.4473, −0.0126 | −1.4473 | no |
| 8.7347 | −1.0031, −0.00018 | −0.00018 | no |
| 10.1330 | 0.0560, 11.588 | 11.588 | **yes: triple crossing** |
| 10.3969 | 0.0999, 0.9317 | 0.9317 | **yes: triple crossing** |

At the two real triple crossings alg = 6. Their Jordan type is not determined. det S⁰ ≠ 0 there (its zeros are at ε ≈ 8.823 and 10.118).

## 8. Verification status and gaps
- Lemma 1 and §5 are SymPy only. A second implementation (e.g. a subresultant computation in the fractions engine, or evaluation of disc and C at the isolating boxes with independent interval code) is planned in ROOT_CONTINUATION_CERTIFICATES (Phase 6).
- Lemma 2 is a finite certificate and the easiest item to re-check independently.
