# RESULT: Jordan type at the real triple crossings, certified by the GEN43 engine

| Field | Value |
|---|---|
| Outcome | **POSITIVE, pointwise.** At each of the four real triple crossings on the witness line, the light eigenvalue of P restricted to ker L has Jordan type exactly **J₃ ⊕ 3J₁** |
| Status in this repository | **CERTIFIED EXTERNALLY**, by the GEN43 engine r10.2.1. That engine is **not** part of this repository. The certificate data are included in `data/gen43_r10_2_1_triple_crossings/`, but this repository contains no script that re-derives them |
| Review status | AI-checked only. Not peer reviewed, no human expert review, not formally verified, not independently replicated |
| Supersedes | The EXPLORATORY entry "real triple crossings at ε ≈ 10.133 and 10.397" in `audit/README.md`, as far as the Jordan type is concerned. The earlier entry is kept unchanged |

## Statement

**Family.** q = (1 + ε)·q_slice, with the parent m, λ, b₁, b₂ and π′ = 1/2. Covector on ray 1, ξ = (s, 0, 1, 0). Write v = ε². The variable η = 35736(q − q_slice), so η² = 2450894·v.

**Where the triple crossings are.**
- Res_u(disc_X Q, C) = κ·v²·p₅(v)², where p₅ is squarefree and of degree 5.
- p₅ has five real roots, in (32, 33), (54, 55), (76, 77), (102, 103) and (108, 109).
- The root near 32.29 is the exact exclusion at ε ≈ 5.68, already recorded in `audit/RESULTANT_AND_ROOT_ISOLATION.md`.
- The roots v ≈ 102.677 and v ≈ 108.096 (ε ≈ 10.133 and 10.397) have u* > 0. These are **real triple crossings**, with algebraic multiplicity 6 on ker L.
- Each comes with both signs of η, which gives **four points**.

**Jordan type at each of the four points: J₃ ⊕ 3J₁.** The ingredients:
1. Algebraic multiplicity is 6. This comes from the ray dispersion identity, the quotient factor and the model identities. It uses the same dispersion factor Q as `data/m2_dispersion.json` (sha256 `baed27fd…`).
2. The geometric multiplicity on ker L is exactly 4:
   - rank [L; P − r] ≥ 12 is certified exactly, by a prime witness in the field K below;
   - ≤ 12 comes from the upper rank identity on ray 1.
3. rank [L; (P − r)²] ≥ 11 is certified exactly by a prime witness, and separately by rigorous intervals. So dim ker (P − r)² on ker L is at most 5.
4. With multiplicity 6 and geometric multiplicity 4, the only possible types are J₃ ⊕ 3J₁ (dim ker (P − r)² = 5) and 2J₂ ⊕ 2J₁ (dim ker (P − r)² = 6). Item 3 rules out 2J₂ ⊕ 2J₁. The engine's 2J₂ route does not pass, as expected.
5. **The number field.** K = ℚ(i)[t, η, s, r]/(p₅(t), η² − 2450894t, s² − u*(t), r² − 1 − u*(t)) is a field of degree 80. Its 80 embeddings cover both triple crossings and both signs of η. Exact rank bounds hold at every embedding.
6. **Constraint solvable.** d₀(η) = det S⁰ ≠ 0 at each point, by rigorous interval enclosure (width about 1e-75).

**Scope.** Four single points (pointwise), with π′ = 1/2. Nothing about neighbourhoods, other parameters, or physics.

## Provenance and checks

| Item | Value |
|---|---|
| Engine | GEN43 r10.2.1, `GEN43_engine_r10_2_1_2026-11-23_1957.zip` (sha256 `56f0a3a8cce07dd5848672f287947f486f09b94ba74ec2c6a95ec2d52e36b5b2`). Not published here |
| Run | the owner's own PC (Ryzen 9 9900X), 2026-11-23, profile r10: 157 results, 0 mathematical failures, 0 software errors |
| Anchor | SHA-256 of the full results-folder manifest: `717c130f94bf8174c39886d3ddfac0a99c6340d2c4d10db59457de63b9c1524a`. The owner stored it outside the results folder and recomputed it on 2026-11-25: MATCH |
| Certificate files | the 34 files in `data/gen43_r10_2_1_triple_crossings/`, byte-identical to the run (re-checked on 2026-11-25 against the run's `certificate_index.json`; subset in `CERTIFICATE_INDEX_SUBSET.json`) |
| Stricter re-check | GEN43 r10.3 "Gate 1b" (2026-11-25), strict gate v2: all four classifications remain VERIFIED_EXACT, with 0 label changes. The 12 scope-migration entries these classifications use were approved by the owner (Sam DeRenzis) on 2026-11-25. That approval is a scope judgement on the software records, not a review of the mathematics |
| Shared inputs | GEN43 transcribes the defining equations (S^ρ, X_A, M, B(ξ), L(ξ), connection and gauge data) from the parent paper and code, and uses a different gamma-matrix representation and its own arithmetic. Agreement with this repository is therefore *shared formulas, different implementation*. It is not independent confirmation |

## Also certified by GEN43 r10.2.1 (summary only; no files included)
The light-factor chain, labelled VERIFIED_EXACT, for every parameter point with det S⁰ ≠ 0, every π′ and every ξ with ξ·ξ ≠ 0:
- L·Y = 0;
- L·P₀ = C₁·L;
- the constrained light factor.

In this repository, LP = C₁L is reproduced by grid proof at π′ = 1/2 (`audit/RESTRICTED_SYMBOL_COMPLETE_PROOF.md`). The GEN43 statement is wider. Its certificates are not included here.

## What this does not establish
- Anything at other points, or on an open set.
- Uniformity, or which part of well-posedness fails.
- Physics.
- Independent or human verification.

A size-3 Jordan block is a stronger defect than the J₂ found on the Milestone 2 family. What it implies for the analytic argument has not been worked out.
