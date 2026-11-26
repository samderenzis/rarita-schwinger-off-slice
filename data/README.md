# data

Exact witnesses and certificates.

| File | What it is | Reader / checker |
|---|---|---|
| `point_A_chain.json`, `point_Am_chain.json`, `point_B_chain.json`, `point_Bm_chain.json` | Off-slice points A, A−, B, B− (see `scripts/points.py`): parameters, covector ξ, eigenvalue r, outcome, ranks, and the exact chain v → w. Entries are `[a, b, c, d]` meaning (a + b·i + (c + d·i)·√δ)/DEN, with δ given in the file. | `scripts/verify_point_fractions.py` (plain Python, fractions only) |
| `m2_minor.json` | Milestone 2 certificate: rows and columns of the minor, the basis K (entries a + b·s), the exact coefficients of d(a, b), and the box certificates. | `scripts/verify_m2_fractions.py` |
| `m2_extended_certificate.json` | Extended Milestone 2 certificate: the exact polynomial matrix M(τ, η), the basis K, and 108 boxes (segment, box, rows, columns, rational C). | `scripts/verify_m2_extended_fractions.py` |
| `m2_gap_certificate.json` | Endpoint-gap certificate: M1 (ray 1), K1, M2 (ray 2), K′, end ε values and 22 boxes. | `scripts/verify_m2_gaps_fractions.py` |
| `m2_dispersion.json` | The ray dispersion factor Q(x, s, ε) (exact coefficients). | `scripts/verify_m2_dispersion_fractions.py` |
| `gen43_r10_2_1_triple_crossings/` | Certificates from the external GEN43 engine r10.2.1 (owner run, 2026-11-23) for the Jordan type J₃ ⊕ 3J₁ at the four real triple crossings: root certificate, exact multimodular ranks in the degree-80 field, interval enclosures, diagnostics and dependency graphs. `CERTIFICATE_INDEX_SUBSET.json` lists each file's SHA-256 from the run. | **None in this repository** (the GEN43 engine is not published here); see `results/RESULT_triple_crossings_GEN43_2026-11-23.md` |
