# rarita-schwinger-off-slice

![Status](https://img.shields.io/badge/status-one--parameter_family_result_·_AI--checked_only-yellow?style=flat-square)
![Scope](https://img.shields.io/badge/scope-exploratory-orange?style=flat-square)

> **This repository is exploratory. Nothing in it is part of the published theorem, and nothing here is established until it is stated as a result, with its evidence and its review status.**

## What this is
This is a follow-up to [rarita-schwinger-well-posedness](https://github.com/samderenzis/rarita-schwinger-well-posedness), the release [v1.0.11](https://github.com/samderenzis/rarita-schwinger-well-posedness/releases/tag/v1.0.11) and [paper version 2](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/main/paper/main.pdf).

That paper proves a local obstruction to loss-free L² well-posedness for a charged, Pauli-coupled massive Rarita–Schwinger field. The proof is at one exact parameter point on the **Pauli slice**, where the charge is q = λ√(2/3)·m.

This repository asks the next question in that project's [ROADMAP](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/main/ROADMAP.md#2-is-the-witness-isolated-or-part-of-a-family), §2:

> **Off the Pauli slice (q ≠ λ√(2/3)·m), near the witness: do nearby characteristic crossings exist, and does the constrained symbol (P restricted to ker L) have a nontrivial Jordan chain at any of them?**

The working plan is in [PLAN.md](PLAN.md). Before any work, read [CONVENTIONS.md](CONVENTIONS.md), which covers sign conventions and the traps inherited from the parent code. Outcomes are reported with [RESULT_TEMPLATE.md](RESULT_TEMPLATE.md). Terms are explained in [GLOSSARY.md](GLOSSARY.md).

## Status
**In short (all AI-checked only; nothing here is released, peer reviewed or human verified):**
- **Established in this repository** (exact certificates; most parts re-checked in a second engine, and the exceptions are listed in [audit/README.md](audit/README.md)):
  - Jordan chains in the constraint kernel at four exact off-slice points.
  - Exactly J₂ ⊕ 3J₁ at both split crossings, for every ε on their full existence ranges. Only q varies; m, λ, b₁, b₂ are fixed.
- **Certified externally:** J₃ ⊕ 3J₁ at the four real triple crossings, by the separate GEN43 engine.
- **Reproduced:** all scripts re-run on the owner's computer, ALL MATCH.
- **Not shown:** anything on an open set of parameters, uniformity in ε, which part of well-posedness fails, or physics. See [What is still open](#what-is-still-open).

**Details, in the order the work was done:**

**PLAN steps 1 and 2.** Step 1: the parent's symbol code now takes q as an independent parameter, and on the slice it reproduces the parent's golden outputs exactly ([log](results/step1_slice_regression_2026-11-22.txt)). Step 2: the reduction identities hold with q free, and at the witness m, λ, b the elimination matrix S⁰ stays invertible for |ε| < 8.82, where q = (1 + ε)·λ√(2/3)·m ([log](results/step2_reduction_q_free_2026-11-22.txt)).

**Four pointwise POSITIVE results (not yet released, AI-checked only).** Off the Pauli slice the parent crossing splits into two crossings. At four exact points, taking both crossings with the charge 0.03% or 0.17% above or below the slice, there is an exact Jordan chain inside the constraint kernel, and the four hypotheses of the parent argument were re-checked at each: [results/RESULT_offslice_points_2026-11-22.md](results/RESULT_offslice_points_2026-11-22.md). These are pointwise witnesses.

**Milestone 2, one direction (not yet released, AI-checked only):** for **every** charge q = (1 + ε)·q_slice with 0 < |ε| ≤ 1, keeping the parent's m, λ, b, both split crossings have a Jordan chain in the constraint kernel. The other hypotheses were re-checked for all such q. The proof is an exact box certificate, checked in two arithmetic engines: [results/RESULT_m2_q_direction_2026-11-22.md](results/RESULT_m2_q_direction_2026-11-22.md). This is a one-parameter family, not an open set in parameter space. **Extended:** the certificate now covers the lower crossing for 0 < |ε| ≤ 2.66 and the upper for 0 < |ε| ≤ 6.28, which is almost their whole existence range: [results/RESULT_m2_extended_2026-11-22.md](results/RESULT_m2_extended_2026-11-22.md). **Closed (development audit):** the remaining gaps are now certified up to and including the exact thresholds ε_L ≈ 2.6674 and ε_U ≈ 6.2912, and at every such crossing the Jordan form is exactly one J₂ plus three J₁: [audit/JORDAN_STRUCTURE_COMPLETE.md](audit/JORDAN_STRUCTURE_COMPLETE.md), [audit/EPSILON_ENDPOINT_CLOSURE.md](audit/EPSILON_ENDPOINT_CLOSURE.md). How the point was found: [notes/STEP4_EXPLORATION.md](notes/STEP4_EXPLORATION.md).

**Triple crossings (certified externally, pointwise):** at the four real triple crossings (ε ≈ ±10.133, ±10.397), the constrained symbol has a Jordan block of size 3 (exact type J₃ ⊕ 3J₁). This was certified by the separate GEN43 engine, which is not part of this repository. Its certificate files are included, but nothing in this repository re-derives them: [results/RESULT_triple_crossings_GEN43_2026-11-23.md](results/RESULT_triple_crossings_GEN43_2026-11-23.md).

**Owner reproduction (2026-11-25):** every script with a committed log was re-run on the owner's own computer (Windows 11, Python 3.14.7) and reproduced its committed output. The result was ALL MATCH; regenerated JSON differs only in Windows line endings. This shows the code is reproducible, not that the mathematics is correct: [results/OWNER_REPRODUCTION_2026-11-25.md](results/OWNER_REPRODUCTION_2026-11-25.md).

## Follow-up: independent certificates (2026-11-27)
The state of this repository on 2026-11-25 is frozen as tag [`frozen-2026-11-25`](https://github.com/samderenzis/rarita-schwinger-off-slice/tree/frozen-2026-11-25). Its exact zip, `rarita-schwinger-off-slice_github_2026-11-25_2230.zip` (SHA-256 `5d92969c8755c3d799910d561379c193b53b114c5a069790ac4f2c265ba362d4`), is attached to release [B6.6.1](https://github.com/samderenzis/rarita-schwinger-off-slice-certificates/releases/tag/B6.6.1) of [rarita-schwinger-off-slice-certificates](https://github.com/samderenzis/rarita-schwinger-off-slice-certificates). That repository checks the zip by hash and builds on it. Since the snapshot, only this section, CHANGELOG.md and SHA256SUMS.txt have changed here.

The certificates repository re-derives several results here in a separate exact-arithmetic verifier, and adds new ones. It is AI-reviewed (ChatGPT, Perplexity) and **not human-verified**.
- **Re-derived, agreeing:** det S⁰ at the witness; the thresholds ε_L and ε_U (exact certificate); J₂ ⊕ 3J₁ at three split crossings inside those ranges (exact certificates, agreeing with Theorem J at those points); the location of the triple crossings. The J₃ ⊕ 3J₁ type at the triple crossings is **not** re-checked there; it still rests on the GEN43 certificate alone.
- **New:** light-cone crossings also exist between the two zeros of det S⁰ (8.823 < |ε| < 10.118) and beyond 10.118, a range not covered by [audit/JORDAN_STRUCTURE_COMPLETE.md](audit/JORDAN_STRUCTURE_COMPLETE.md). All 20 sampled there show J₂ ⊕ 3J₁ numerically, and the one at ε ≈ 9.958 is certified exactly. These points lie outside C₀, the region around the slice bounded by det S⁰ = 0, where the characteristic roots can be non-real (certified at two points).
- **New:** near the slice, all characteristic roots are real on certified boxes covering about 97.5 % of C₀ for 0.10 ≤ √6·λ·m ≤ 0.17. Real roots give weak hyperbolicity only.
- **Still not shown:** an open set in (m, λ, b, q), uniformity in ε, which part of well-posedness fails, physics, and human review.

## What is still open
Tracked in [PLAN.md](PLAN.md):
- **Varying m, λ, b₁, b₂.** At present every family result holds along one line in parameter space (q varies). Whether the obstruction persists on an open set needs its own proof.
- **Uniformity in ε** of the radius and constants. This is needed for a common ball.
- **How severe the failure is.** Derivative loss or worse: the parent ROADMAP, §3.
- **Independent (non-AI) review** of everything here.

## How to check it
- **Checksums.** From the folder that contains the repository, run `sha256sum -c rarita-schwinger-off-slice/SHA256SUMS.txt`. On Windows, use `certutil -hashfile <file> SHA256` for single files.
- **One script.** From the repository root, run `python3 -I scripts/<script>.py` with SymPy 1.14.0 and mpmath 1.3.0. The `verify_*_fractions.py` second-engine scripts need only the standard library. [scripts/README.md](scripts/README.md) gives the command and the expected run time for each script.
- **Everything.** From the folder that contains the repository, run:
  ```
  python3 -I rarita-schwinger-off-slice/results/owner_reproduction_2026-11-25/offslice_rerun.py --repo rarita-schwinger-off-slice --parent <path to rarita-schwinger-well-posedness v1.0.11>
  ```
  This re-runs every script on a copy and compares each output with the committed one. It takes about 15–35 minutes. See [CONVENTIONS.md](CONVENTIONS.md) for the known traps, such as Windows line endings.

## Ground rules (carried over from the parent project)
- **Exact arithmetic for anything called a result.** High-precision numerics are used only to decide where exact work is worth doing.
- **SHA-256 for every file**, in `SHA256SUMS.txt`.
- **Failed attempts and negative results stay on the record**, in [notes/FAILED_ATTEMPTS.md](notes/FAILED_ATTEMPTS.md).
- **Outcomes are reported separately:** positive, negative, and "no crossing found". A pointwise result is not a family result. Persistence on an open set needs its own proof.
- **No claims of novelty, independent replication, human verification or physical consequence** unless they actually exist.
- This is AI-collaborative work. Its review status is stated plainly.

## Layout
| Folder | What goes there |
|---|---|
| [`scripts/`](scripts/) | computation scripts, with q as an independent parameter (see [scripts/README.md](scripts/README.md)) |
| [`data/`](data/) | exact witnesses and certificates, including the externally certified triple-crossing certificates (see [data/README.md](data/README.md)) |
| [`audit/`](audit/) | written proofs and audit documents (coercivity, restricted symbol, Jordan structure, endpoint closure) with their scripts; claim-status changes are listed in [audit/README.md](audit/README.md) |
| [`results/`](results/) | reference outputs, and the logs of runs |
| [`notes/`](notes/) | working notes, and failed attempts |

## License
The code is under the [MIT License](LICENSE). Written material (notes, and any future paper) is under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact
Sam DeRenzis · waitandhope123@gmail.com · Comments and corrections are welcome: please open an issue.
