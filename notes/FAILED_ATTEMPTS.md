# Failed attempts and negative results

Every attempt that did not work is kept here: a wrong approach, a numerical dead end, a crossing that turned out semisimple, a search that found nothing. Each entry says:
- what was tried;
- what happened;
- what it does **not** rule out.

Nothing is deleted. Entries are in the order they happened (F-1, F-2, …); new entries go at the end.

## F-1 (2026-11-22) Rational ε with a rational crossing
- **Tried:** rational ε (q a rational multiple of the slice value) such that the crossing u = |ξ_a|²/|ξ_b|² is rational.
- **Happened:** none found. The discriminant of the crossing polynomial in u is 1646·(square)·ε²·(91809 − 823ε²), so a rational crossing needs 1646(91809 − 823ε²) to be a rational square. A quick local argument at the prime 823 (823 ≡ 7 mod 8) suggests this has no rational solution; that argument has not been written up or checked carefully.
- **Does not rule out:** crossings at rational ε (they exist, at quadratic irrational u); rational crossings for irrational ε (point A is one).

## F-2 (2026-11-22) Pythagorean covector on the ray ξ = (s, 0, 1, 0)
- **Tried:** s and |ξ| both rational (s = (1 − τ²)/(2τ), τ rational with denominator < 400), with rational ε².
- **Happened:** only τ = 9/11 (the parent s* = 20/99, with the far root ε ≈ 9.96) was found. The conditions combine to a curve of genus 1, so rational points may be rare.
- **Does not rule out:** such points of larger height. The search moved to general covectors (a₁, a₂, b₁, b₂), which gave point A.

## F-3 (2026-11-22) Conic search returned nothing (bug, fixed)
- **Tried:** the first conic search near u* used a `lambdify`-compiled formula.
- **Happened:** it silently produced floating-point values (denominators 2⁵⁷), so the "is a square" tests never passed and the search reported zero candidates. Found by checking the denominators; replaced by exact substitution, which found point A. The fixed search is part D/E of `scripts/step4_dispersion_and_crossings.py`.
- **Does not rule out:** nothing; this was a code error, kept here as a record.

## F-4 (2026-11-22) Symbol assumptions in step 4 (bug, fixed)
- **Tried:** substituting x = √(1+u) into the characteristic polynomial.
- **Happened:** the script's `x` had `real=True`, while SymPy's `charpoly(x).as_expr()` returns a plain `x`, so the substitution silently did nothing and a later step crashed. Fixed by using a plain symbol. No reported number was affected (the crash happened before any output was used).

## F-5 (2026-11-22) Milestone 2 box centred at the wrong point (bug, fixed)
- **Tried:** the first Milestone 2 box certificate used the parametrisation s = 2τ/(1 − τ²) but centred the box at τ = 9/11. That value belongs to the other parametrisation, s = (1 − τ²)/(2τ), from F-2.
- **Happened:** τ = 9/11 is s ≈ 4.95, far from any crossing. The "crossings stay in the box" test passed vacuously (no crossing was near), and the minor was certified at a meaningless point. Caught because the box half-widths did not grow with ε₀ and the result looked too good. A check then showed the crossing polynomial is nonzero at the claimed centre.
- **Fix:** centre at τ = 1/10 (s = 20/99 exactly). The script now asserts that the centre maps to s* = 20/99 and r* = 101/99, and that the crossing polynomial vanishes there. Nothing from the wrong run was reported.
- **Does not rule out:** nothing; a code error, kept as a record.

## F-6 (2026-11-22) Second check of the extended certificate with naive interval arithmetic
- **Tried:** re-checking the 108 rank-8 boxes with rational interval arithmetic on each monomial (a different bound from the planner's Taylor bound), with up to 4 levels of subdivision.
- **Happened:** 62 boxes confirmed, 46 not. The η ranges are wide (η is up to about 10⁴), and naive interval powers lose the cancellation that the centred Taylor form keeps. A bound that is too crude is not evidence against the claim.
- **Replaced by:** the same centred Taylor bound in an independent implementation (Taylor shift by synthetic division instead of binomial expansion), which confirms all 108 boxes. Negative controls fail as they should.
- **Does not rule out:** nothing; the weaker bound is kept as a record.

## F-7 (2026-11-22) Second check of the gap certificate: wrong zero test (bug in the checker, fixed)
- **Tried:** the first version of `scripts/verify_m2_gaps_fractions.py` tested "one nonzero entry per row of K' outside block 2" by counting non-empty coefficient lists.
- **Happened:** SymPy records a zero polynomial as a single zero term, so every entry counted as nonzero and the check FAILED although K' was correct. Found by printing the stored entries.
- **Fix:** test the coefficients themselves; also check that each such entry is exactly (1 + t²) (that is, K' = 1 there) and that the positions form the 12 × 12 identity. The certificate itself was unchanged.
- **Does not rule out:** nothing; a checker bug, kept as a record. It shows the checker does fail loudly on a malformed condition.

## F-8 (2026-11-25) Owner reproduction: first runner reported false differences (bug in the runner, fixed)
- **Tried:** the first version (1.0) of the reproduction runner compared each regenerated data file with the committed one byte for byte (SHA-256). The run was on the owner's Windows 11 machine with Python 3.14.7.
- **Happened:** 5 jobs were reported DIFFER. They were `step5_point_exact.py` for A, A−, B and B−, and `m2_family_q_direction.py`, which are exactly the scripts that write indented JSON in text mode. On Windows those files get CRLF line endings, so their bytes and SHA-256 change while the content does not. The 3 scripts that write single-line JSON matched byte for byte. The second-engine checks of the committed data all matched. Step 1 was NOT RUN in that attempt because the parent repository was not given.
- **Fix:** runner 1.1 compares the JSON after converting CRLF to LF. It reports such a job as MATCH-EOL, maps the hash printed in the log only when the converted file equals the committed one, and still reports any real content change as DIFFER. Before use it was tested by converting files to CRLF (MATCH-EOL) and by CRLF plus one changed character (DIFFER). The re-run with 1.1 and the parent repository gave ALL MATCH (`results/OWNER_REPRODUCTION_2026-11-25.md`). No script, result or data file was changed.
- **Does not rule out:** nothing about the mathematics. It records that the JSON outputs are not byte-portable across operating systems (CONVENTIONS.md, "JSON line endings").
