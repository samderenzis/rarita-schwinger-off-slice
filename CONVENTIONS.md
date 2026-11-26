# Conventions

These are the conventions and known traps carried over from the parent repository, [rarita-schwinger-well-posedness](https://github.com/samderenzis/rarita-schwinger-well-posedness), release v1.0.11. File and line references below are to that release. (Earlier versions of this file cited v1.0.9. The three files referred to, `scripts/bhmap_defs.py`, `scripts/witness_min.py` and `bhmap_validate_v1_7.py`, are byte-identical in v1.0.9 and v1.0.11 (same SHA-256), and the quoted lines were re-checked in v1.0.11 on 2026-11-22.)

## Signs and names
- **The matrix M.**
  - In the paper, M is Γ^{0ak}.
  - In the parent code, the variable `Mmat` is the **negative** of the paper's M. The paper's M is `MMmat = -Mmat`: see [`scripts/bhmap_defs.py`](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/v1.0.11/scripts/bhmap_defs.py), line 54, and the validator, line 120, both commented "MM^{ak} = Gamma^{0ak}".
  - Any formula copied between the paper and code here must account for this sign.
- **Parameters.** These are the names used in the parent code:

  | Name | Meaning |
  |---|---|
  | λ (`lam`) | Pauli coupling |
  | m | mass |
  | q | charge |
  | b₁, b₂ | magnetic field components |
  | π′ | constraint parameter |
  | ξ | covector |
  | r | eigenvalue (characteristic speed) |

## The Pauli slice is hard-wired in the parent code
- [`scripts/witness_min.py`](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/v1.0.11/scripts/witness_min.py), line 12, sets `q = lam*sqrt(2/3)*kap*m`.
- [`scripts/bhmap_defs.py`](https://github.com/samderenzis/rarita-schwinger-well-posedness/blob/v1.0.11/scripts/bhmap_defs.py), line 33, sets `q = lam*sqrt(2/3)*m`, commented "Pauli slice".
- **The factor `kap` does not leave the slice.** In witness_min.py, `kap` scales λ in both q and the Pauli term (line 16), so q = λ√(2/3)·m still holds.
- Code in this repository must take **q as an independent parameter**. Before such code is used off the slice, it must reproduce the parent's golden outputs exactly when q is set back to λ√(2/3)·m.

## Arithmetic
- **Results use exact arithmetic only.**
  - The parent witness lives over ℚ(i).
  - Off the slice, the parent's ROADMAP §2 expects crossings at nested-radical points. Ranks and chains must then be computed over the corresponding algebraic number field, stated explicitly.
- **High-precision floating point is for exploration only**: deciding where exact work is worth doing. A floating-point computation is never reported as a result.
- Each exact claim states:
  - the field it was computed over;
  - the engine or engines used;
  - whether it was cross-checked in a second engine.

## Files and records
- Files are listed in `SHA256SUMS.txt` (format: `hash  path`, sorted). The file does not list itself. Its paths start with `rarita-schwinger-off-slice/`, so check it from the folder that **contains** the repository: `sha256sum -c rarita-schwinger-off-slice/SHA256SUMS.txt`. (Earlier versions of this file called it `SHA256SUMS`.)
- Runs are logged in `results/` with:
  - the date;
  - the Python and package versions;
  - the command;
  - the exit code;
  - the hashes of the inputs and outputs.
- Failed attempts go in [`notes/FAILED_ATTEMPTS.md`](notes/FAILED_ATTEMPTS.md). Nothing is deleted. A correction is a new entry that refers to the old one.
- Outcomes are reported using [RESULT_TEMPLATE.md](RESULT_TEMPLATE.md).

## Wording
- A pointwise result is called a pointwise result.
- "Family", "open set" and "persists" are used only with a proof that supports them.
- Agreement between this project's own scripts is code consistency, not independent confirmation.
- Principal-symbol characteristic speeds are not signal speeds.

## Off-slice parametrisation
- q = (1 + ε)·q_slice with q_slice = λ√(2/3)·m. Unless a file says otherwise, m, λ, b₁ and b₂ keep the parent witness values. The scripts also use η = 35736·(q − q_slice), so η² = 2450894·ε².
- Crossings are tracked along **ray 1**, ξ = (s, 0, 1, 0), and **ray 2**, ξ′ = (1, 0, t, 0). The other terms (light eigenvalue, split and triple crossings, ε_L, ε_U) are defined in [GLOSSARY.md](GLOSSARY.md).

## Running the scripts
- Run every script from the repository root with `python3 -I`. Use SymPy 1.14.0 and mpmath 1.3.0, the versions the logs were made with. The `verify_*_fractions.py` scripts need only the Python standard library.
- **Logs record script hashes.** Every log in `results/` records the SHA-256 of the script that produced it, and of helper scripts it imports. Changing a script, even only its formatting, makes its committed logs out of date. A script change must therefore come with regenerated logs and a CHANGELOG entry. Until then, the old logs stay as they are and the difference is noted.
- **JSON line endings.** The scripts that write indented JSON (`step5_point_exact.py`, `m2_family_q_direction.py`) open their output in text mode. On Windows the output gets CRLF line endings, so it is byte-different from the committed LF files although the content is identical. Compare after converting CRLF to LF. `results/owner_reproduction_2026-11-25/offslice_rerun.py` does this. The scripts were left unchanged for the reason in the previous point.
- **The `m2_family_extended` log predates the current `m2_box_tools.py`.** The log records an earlier hash of the helper. The certificate it writes is byte-identical with the current helper (`results/OWNER_REPRODUCTION_2026-11-25.md`).

## Results certified outside this repository
- A result certified by an engine that is not in this repository is labelled **CERTIFIED EXTERNALLY**. Its certificate files are included unchanged, together with their recorded SHA-256 values and the provenance of the run (engine release and hash, who ran it and when, and how the files were checked).
- Such a result is never described as re-derivable from this repository. It is not combined with the repository's own results into a stronger claim. It carries the engine's shared-input caveats, for example a common transcription of the defining equations.
- Current case: the triple-crossing Jordan type, `results/RESULT_triple_crossings_GEN43_2026-11-23.md`.
