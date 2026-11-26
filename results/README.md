# results

Reference outputs of the scripts and the logs of runs. Each run records the date, the tool versions, the command, the exit status and the hashes of its inputs.

| File | What it is |
|---|---|
| `RESULT_triple_crossings_GEN43_2026-11-23.md` | **POSITIVE, pointwise, certified externally:** J₃ ⊕ 3J₁ at the four real triple crossings, by GEN43 r10.2.1 (engine not in this repository; certificates in `data/gen43_r10_2_1_triple_crossings/`). |
| `OWNER_REPRODUCTION_2026-11-25.md`, `owner_reproduction_2026-11-25/` | Owner re-run of every script with a committed log (Windows 11, Python 3.14.7): ALL MATCH. Reproducibility only, not correctness. |
| `step1_slice_regression_2026-11-22.txt` | PLAN step 1: the q-independent code reproduces the parent (v1.0.11) exactly on the Pauli slice. Exit code 0 (ALL PASS). Regression only, not a result. |
| `m2_endpoint_gaps_2026-11-22.txt`, `verify_m2_gaps_fractions_2026-11-22.txt` | Endpoint-gap certificate and its second check (exit 0); see `audit/EPSILON_ENDPOINT_CLOSURE.md`. |
| `audit_*_2026-11-22.txt` | Logs of the `audit/` scripts (coercivity constants, restricted symbol, Jordan classification, endpoint multiplicity), all ALL PASS. |
| `regression_verify_m2_extended_new_boxtools_2026-11-22.txt` | The extended certificate re-verified with the modified `m2_box_tools.py` (ALL PASS). |
| `RESULT_m2_extended_2026-11-22.md` | **POSITIVE, extended family:** Jordan chains at the lower crossing for 0 < \|ε\| ≤ 2.66 and at the upper for 0 < \|ε\| ≤ 6.28 (108 certified segments). |
| `m2_family_extended_2026-11-22.txt`, `verify_m2_extended_fractions_2026-11-22.txt` | Run logs of the extended certificate and its second check (exit 0). |
| `m2_dispersion_export_2026-11-22.txt`, `verify_m2_dispersion_fractions_2026-11-22.txt` | Ray dispersion relation: SymPy export of Q and the fractions-only proof of the full characteristic polynomial (exit 0). |
| `RESULT_m2_q_direction_2026-11-22.md` | **POSITIVE, one-parameter family:** Jordan chains at both crossings for every 0 < \|ε\| ≤ 1 (q varies; m, λ, b fixed), with all hypotheses re-checked. Uses RESULT_TEMPLATE. |
| `m2_rotation_invariance_2026-11-22.txt`, `m2_family_q_direction_2026-11-22.txt`, `verify_m2_fractions_2026-11-22.txt`, `m2_curved_symbol_general_2026-11-22.txt` | Run logs of the Milestone 2 scripts (all exit 0). |
| `RESULT_offslice_points_2026-11-22.md` | **POSITIVE (pointwise) at four points** A, A−, B, B−: exact constrained Jordan chains at off-slice crossings, with all four hypotheses re-checked. Uses RESULT_TEMPLATE. |
| `step4_dispersion_and_crossings_2026-11-22.txt` | Steps 3–4: exact dispersion relation along a ray with ε symbolic, the crossing discriminant, a numerical Jordan probe (exploration only), and the exact search that produced points A and B. |
| `step5_point_<P>_exact_2026-11-22.txt` (P = A, Am, B, Bm) | Step 5 at each point, SymPy over ℚ(i, √δ): det S⁰, charpoly, ranks 4, 5, 5, chain, π′-independence. Exit 0. |
| `verify_point_<P>_fractions_2026-11-22.txt` | Second engine (pure-Python fractions) at each point; also reproduces the parent chain on the slice. Exit 0. |
| `step5_point_<P>_curved_symbol_2026-11-22.txt` | Curved-symbol identification at each point (q symbolic). Exit 0. |
| `step2_reduction_q_free_2026-11-22.txt` | PLAN step 2: reduction identities (divergence identity, constraint intertwining, π′-factorisation, invariance of ker L) with q free; det S⁰ as an exact function of q. Exit code 0 (ALL PASS). Algebraic identities, SymPy only. |
