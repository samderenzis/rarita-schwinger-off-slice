# scripts

Requirements: Python 3, SymPy 1.14.0 (the version the parent project's golden outputs were produced with). Run every script with `python3 -I`.

| Script | What it does |
|---|---|
| `rs_model.py` | The parent's `bhmap_defs.py` (v1.0.11) with the charge **q as an independent parameter**: `Model(m, lam, b1, b2, pip, q=None)`. `q=None` gives the Pauli slice q = λ√(2/3)·m, exactly as in the parent. Nothing else is changed; the docstring lists the differences. |
| `parent_shim_runner.py` | Runs one **unmodified** parent check script with `rs_model` standing in for `bhmap_defs`, passing q explicitly as λ√(2/3)·m. Used by the regression test. |
| `witness_check.py` | The parent's `witness_min.py` checks, rebuilt on `rs_model` with q passed explicitly. |
| `step2_reduction_q_free.py` | PLAN step 2 test. With q an independent symbol: divergence identity and constraint intertwining on the curved background (arbitrary fields, symbolic parameters); π′-factorisation and LP = C₁L at the witness m, λ, b with symbolic ξ, π′, q; det S⁰ as an exact function of q and of ε. About 60 seconds. Run: `python3 -I scripts/step2_reduction_q_free.py --out results/step2_reduction_q_free_<date>.txt` |
| `step4_dispersion_and_crossings.py` | Steps 3–4: q = (1+ε)·q_slice. Exact dispersion relation on ker L along ξ = (s, 0, 1, 0) with ε symbolic, the crossing polynomial and its discriminant, 60-digit numerical Jordan probe at the crossings (exploration only), and the exact search that finds point A. ~20 s. |
| `points.py` | Exact constants of the off-slice points A, A−, B, B− (fractions only), how they were chosen, and a self-check (`python3 -I scripts/points.py`). |
| `step5_point_exact.py` | Step 5 at one point (`--point A`, `A-`, `B`, `B-`), over ℚ(i, √δ) (SymPy DomainMatrix): det S⁰, characteristic polynomial, ranks of (P − r)^j on ker L, the outcome (positive or negative, not assumed), an explicit chain if positive, π′-independence; writes `data/point_<P>_chain.json`. ~2 min. |
| `verify_point_fractions.py` | Second engine: pure Python fractions with its own arithmetic in ℚ(i, √δ), no SymPy. Rebuilds Γ, S^ρ, M, L, re-checks the chain and the ranks against the recorded outcome, and reproduces the parent chain on the slice. ~3 s. |
| `step5_curved_symbol.py` | Curved-symbol identification at one point (`--point`), parent check_5 analogue, with q symbolic. ~12 s. |
| `m2_rotation_invariance.py` | Milestone 2: exact rotation equivariance of P and L within each sphere factor (symbolic ξ, q, angle). ~10 s. |
| `m2_family_q_direction.py` | Milestone 2: crossings real for 0 < ε² ≤ 1, the 8×8 minor d(a, b) of Dn·(P − \|ξ\|)K, exact Taylor-box certificates and box containment, so Jordan chains at both crossings for every 0 < \|ε\| ≤ 1; writes `data/m2_minor.json`. ~30 s. |
| `verify_m2_fractions.py` | Second engine for the Milestone 2 certificate (pure fractions): 1089-point polynomial identity of the minor, box inequalities, Sturm containment. ~2 min. |
| `m2_box_tools.py` | Exact (fractions) tools for the extended certificate: Taylor-box rank-8 bound ‖CA − I‖_∞ < 1, crossing polynomial in τ and in v = η², Sturm counts, segment containment with branch identification. |
| `m2_family_extended.py` | Extended Milestone 2 certificate: exact thresholds ε_L, ε_U; 108 ε-segments covering the lower crossing to \|ε\| = 2.66 and the upper to 6.28, both signs; writes `data/m2_extended_certificate.json`. ~20 s. |
| `verify_m2_extended_fractions.py` | Second check without SymPy: M proved on a 5×5 grid from the independent construction of P, rank bounds re-done with an independent Taylor-shift implementation, containment and coverage re-checked. ~10 s. |
| `m2_endpoint_gaps.py` | Closes the endpoint gaps: 22 boxes (ray 1 lower to 1067/400 ≥ ε_L; ray 2 upper to 62913/10000 ≥ ε_U), end boxes touching the endpoint direction. Writes `data/m2_gap_certificate.json`. ~11 s. |
| `verify_m2_gaps_fractions.py` | Second check of the gap certificate without SymPy (M1, M2 proved on grids; K′ checks; independent Taylor shift). ~5 s. |
| `m2_dispersion_export.py` | Exports the exact ray dispersion factor Q (P on ker L, ε symbolic) to `data/m2_dispersion.json`. ~9 s. |
| `verify_m2_dispersion_fractions.py` | Fractions-only proof of the full 16×16 characteristic polynomial along the ray (17 × 65 grid with degree bounds). ~12 min on 2 cores (`NPROC` sets the number of processes). |
| `m2_curved_symbol_general.py` | Curved-symbol identification for every covector and every q (symbolic). ~8 s. |
| `step1_slice_regression.py` | PLAN step 1 test. Checks the parent files' SHA-256, then shows that on the slice the new code reproduces the parent's golden outputs line for line, and the parent's S^ρ, X, M, P(ξ\*), L(ξ\*) entry by entry. Also checks that q is really free (S^ρ is affine in q and moves when q moves). Exit code 0 only if everything passes. |

Run the step 1 test against a copy of the parent repository at v1.0.11:

```
python3 -I scripts/step1_slice_regression.py --parent <path to rarita-schwinger-well-posedness> --out results/step1_slice_regression_<date>.txt
```

It takes about 70 seconds; most of it is the parent's `check_4_constraint_intertwining.py`.

Passing this test means the new code agrees with the parent code on the slice. That is code consistency, not independent confirmation, and it says nothing about any point off the slice.
