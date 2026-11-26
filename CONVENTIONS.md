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
- Released files are listed in `SHA256SUMS` (format: `hash  path`, sorted).
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
