# Owner reproduction, 2026-11-25

Every script in this repository that has a committed log was re-run on the owner's own computer. Each new output was compared with the committed output.

**Result: ALL MATCH (exit code 0).**

## Run details

| Field | Value |
|---|---|
| Machine | owner's PC: Windows 11, Python 3.14.7, SymPy 1.14.0, mpmath 1.3.0 |
| Original logs | produced on Linux with Python 3.13.18 |
| Repository state | the files of this commit before this record was added. `SHA256SUMS.txt` had sha256 `8aa388b53cea9f057b687d33bce65fc2c3724772d2c1c28000e13608ece81ef9`, and 89/89 files were OK |
| Parent | rarita-schwinger-well-posedness v1.0.11 (step 1; the step 1 script pins the parent file hashes itself) |
| Runner | `owner_reproduction_2026-11-25/offslice_rerun.py` (offslice-owner-rerun 1.1). It works on a copy and changes nothing |
| Records | `owner_reproduction_2026-11-25/SUMMARY.txt` and `MANIFEST.json` (per-job status, log hashes, data hashes) |

## Outcome (31 jobs)

- **25 MATCH.** The new log is identical to the committed log. Only these are ignored: date, Python version, command line, timing, and file paths of outputs.
- **5 MATCH-EOL.** These jobs write indented JSON: `step5_point_exact` for A, A−, B and B−, and `m2_family_q_direction`. On Windows, Python writes those files with CRLF line endings. After converting CRLF to LF, each regenerated file is byte-identical to the committed one in `data/`.
- **1 MATCH-KNOWN.** `results/m2_family_extended_2026-11-22.txt` records `scripts/m2_box_tools.py` with sha256 `b5ac2af0…`. The current file is `4c6ed9b6…`, because the tools were extended for the endpoint gaps (see CHANGELOG). The certificate `data/m2_extended_certificate.json` regenerated with the current tools is byte-identical (sha256 `7bcf43c6…`). `results/regression_verify_m2_extended_new_boxtools_2026-11-22.txt` already records that the certificate re-verifies with the new tools.

Total time: 879 seconds.

An earlier attempt the same evening (runner 1.0, without the parent) reported the 5 MATCH-EOL jobs as DIFFER. Its byte comparison did not allow for Windows line endings. Runner 1.1 fixes that, and still flags any real content change.

## What this shows, and what it does not

- **It shows** that the scripts are deterministic and reproduce their committed outputs on a second machine, operating system and Python version.
- **It is not** a proof of mathematical correctness.
- **It is not** independent review or human verification. The results in this repository remain **AI-checked only**.

## Known, not yet changed
- The JSON outputs are not written with `newline='\n'`, so on Windows they come out byte-different (CRLF). The scripts were left unchanged, because every committed log records the SHA-256 of the script that produced it.
- The `m2_family_extended` log predates the current `m2_box_tools.py`; see above.
