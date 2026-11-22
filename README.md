# rarita-schwinger-off-slice

![Status](https://img.shields.io/badge/status-planning_·_no_results_yet-lightgrey?style=flat-square)
![Scope](https://img.shields.io/badge/scope-exploratory-orange?style=flat-square)

> **This repository is exploratory. Nothing in it is part of the published theorem, and nothing here is established until it is stated as a result, with its evidence and its review status.**

## What this is
This is a follow-up to [rarita-schwinger-well-posedness](https://github.com/samderenzis/rarita-schwinger-well-posedness), the release [v1.0.9](https://github.com/samderenzis/rarita-schwinger-well-posedness/releases/tag/v1.0.9) and [paper version 2](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/main/paper/main.pdf).

That paper proves a local obstruction to loss-free L² well-posedness for a charged, Pauli-coupled massive Rarita–Schwinger field. The proof is at one exact parameter point on the **Pauli slice**, where the charge is q = λ√(2/3)·m.

This repository asks the next question in that project's [ROADMAP](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/main/ROADMAP.md#2-is-the-witness-isolated-or-part-of-a-family), §2:

> **Off the Pauli slice (q ≠ λ√(2/3)·m), near the witness: do nearby characteristic crossings exist, and does the constrained symbol (P restricted to ker L) have a nontrivial Jordan chain at any of them?**

The working plan is in [PLAN.md](PLAN.md).

## Status
Planning. No computations have been run here yet.

## Ground rules (carried over from the parent project)
- **Exact arithmetic for anything called a result.** High-precision numerics are used only to decide where exact work is worth doing.
- **SHA-256 for every released file**, in `SHA256SUMS`.
- **Failed attempts and negative results stay on the record**, in [notes/FAILED_ATTEMPTS.md](notes/FAILED_ATTEMPTS.md).
- **Outcomes are reported separately:** positive, negative, and "no crossing found". A pointwise result is not a family result. Persistence on an open set needs its own proof.
- **No claims of novelty, independent replication, human verification or physical consequence** unless they actually exist.
- This is AI-collaborative work. Its review status is stated plainly.

## Layout
| Folder | What goes there |
|---|---|
| [`scripts/`](scripts/) | computation scripts, with q as an independent parameter |
| [`data/`](data/) | exact witnesses and certificates, if any are found |
| [`results/`](results/) | reference outputs, and the logs of runs |
| [`notes/`](notes/) | working notes, and failed attempts |

## License
The code is under the [MIT License](LICENSE). Written material (notes, and any future paper) is under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact
Sam DeRenzis · waitandhope123@gmail.com · Comments and corrections are welcome: please open an issue.
