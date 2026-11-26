# audit

Development-phase audit of the off-slice results (2026-11-22), written in response to the reviewer instructions of that date. Nothing here is released. AI-checked only.

## Documents
| File | Contents |
|---|---|
| `OFFSLICE_COERCIVITY_FINAL.md` | **GEN43 Phase 1 final.** Audit below plus Appendix B constants (N₁, N₂) off the slice and a pointwise / compact-set / full-axis table. |
| `RESTRICTED_SYMBOL_COMPLETE_PROOF.md` | **GEN43 Phase 1 final.** R1–R8 with a summary theorem; rank ≥ 8 coverage now includes ray 2 and the end boxes. |
| `RESULTANT_AND_ROOT_ISOLATION.md` | **GEN43 Phase 1.** Res_u(disc_X Q, C) = κv²p₅², mod-11 irreducibility certificate, isolating intervals, branch identification at v₀ ≈ 32.29. |
| `JORDAN_STRUCTURE_COMPLETE.md` | **GEN43 Phase 1 final.** Theorem J on the full intervals 0 < \|ε\| ≤ ε_L (lower) and ≤ ε_U (upper), endpoints included. |
| `EPSILON_ENDPOINT_CLOSURE.md` | **GEN43 Phase 1.** Ray 2, basis K′, end boxes and the continuity argument, endpoint multiplicity, margins. |
| `OFFSLICE_COERCIVITY_AUDIT.md` (superseded draft, kept) | Every step of the paper's analytic chain checked for use of the slice relation (one use: the parameter choice); exact q-dependence of the coercivity constants; standalone off-slice coercivity proposition with explicit r₀(q) and q-independent c₀; off-slice transfer theorem. |
| `RESTRICTED_SYMBOL_PROOF.md` (superseded draft, kept) | ker L invariance → block form → restricted characteristic polynomial → algebraic and geometric multiplicity, with every denominator and factor accounted for. |
| `JORDAN_STRUCTURE_CLASSIFICATION.md` (superseded draft, kept) | Theorem J: exactly J₂ ⊕ 3J₁ at every certified crossing; where the structure changes. |

## Scripts and logs
| Script | Engine | Log |
|---|---|---|
| `coercivity_constants.py` | SymPy (exact) | `results/audit_coercivity_constants_2026-11-22.txt` |
| `restricted_symbol_fractions.py` | fractions engine (no SymPy) | `results/audit_restricted_symbol_fractions_2026-11-22.txt` |
| `jordan_classification.py` | SymPy (exact) | `results/audit_jordan_classification_2026-11-22.txt` |
| `endpoint_multiplicity.py` | SymPy (exact) | `results/audit_endpoint_multiplicity_2026-11-22.txt` |
| `../scripts/m2_endpoint_gaps.py` | SymPy builds M; exact fractions for all inequalities | `results/m2_endpoint_gaps_2026-11-22.txt` |
| `../scripts/verify_m2_gaps_fractions.py` | fractions engine (no SymPy) | `results/verify_m2_gaps_fractions_2026-11-22.txt` |

Run from the repository root: `python3 -I audit/<script>.py --out results/<log>.txt`.

## Claims whose status changed (running list)
| Claim | Before | Now | Why |
|---|---|---|---|
| Coercivity (Prop. 4.2) holds off the slice | SUPPORTED (proof reading) | ESTABLISHED | written standalone proof; the only q-dependent constant computed exactly: m₀(q) ≤ √(10 + 10λ²/3 + 16m² + 72q²); c₀ q-independent |
| Analytic chain off the slice (Theorem T) | SUPPORTED | ESTABLISHED (witness m, λ, b; det S⁰ ≠ 0) | dependency table; Proposition N |
| Invariance of ker L (LP = C₁L) | ESTABLISHED (SymPy) | REPRODUCED (π′ = 1/2) | fractions-engine grid proof |
| Algebraic multiplicity ≥ 5 at crossings | ESTABLISHED (SymPy) | REPRODUCED | full charpoly grid proof plus quotient factor (R3–R7) |
| "One J₂" at the crossings | observed at 4 exact points and in numerics | ESTABLISHED on the certified ranges (exactly J₂ ⊕ 3J₁) | alg = 5 exactly (resultant + exact exclusion at ε ≈ 5.68), geo = 4 exactly (rank ≤ 8 identically, ≥ 8 by boxes) |
| Light eigenvalue semisimple off the crossings | — | ESTABLISHED (new) | alg = 4 and rank ≤ 8 identically |
| Real triple crossings at ε ≈ 10.133, 10.397 | — | EXPLORATORY (new) | common real zeros of disc_X Q and C, outside the family |
| Jordan type at the real triple crossings (4 points, π′ = 1/2) | EXPLORATORY | J₃ ⊕ 3J₁, CERTIFIED EXTERNALLY (GEN43 r10.2.1; added 2026-11-25) | exact ranks in a degree-80 field; see `results/RESULT_triple_crossings_GEN43_2026-11-23.md` |
| Jordan chain (geo = 4) on the gaps 2.66 < \|ε\| ≤ ε_L (lower), 6.28 < \|ε\| ≤ ε_U (upper) | INCOMPLETE | ESTABLISHED (second check shares containment code) | 22 boxes incl. end boxes containing the endpoint directions; ray 2 with hand-pivoted K′ |
| Algebraic multiplicity 5 at the endpoint crossings ε = ±ε_L, ±ε_U | — | ESTABLISHED (SymPy only, new) | homogenised Q, X = 1 simple root |
| Theorem J range | 0 < \|ε\| ≤ 2.66 / 6.28 | 0 < \|ε\| ≤ ε_L / ε_U (full existence intervals, endpoints included) | the two lines above |
| p₅ irreducible | ESTABLISHED (SymPy factor) | ESTABLISHED with a finite certificate (mod 11) | easy to re-check by hand |
| Appendix B constants off the slice | finite for each q (reading) | ESTABLISHED with explicit q-growth (SymPy only) | ‖∂M₀/∂cᵢ‖ exact |
