# Result: exact constrained Jordan chains at four off-slice crossings (points A, A−, B, B−)   (date: 2026-11-22)

Four pointwise results, each POSITIVE. They share the parameters, the field type, the method and the hypothesis checks, so one report covers them. Each point is a separate pointwise statement.

```
Outcome type: POSITIVE (constrained Jordan chain found) — at each of the four points separately

Parameter points (common: λ² = 1/2, m² = 823/142944, b₁ = 1, b₂ = 1/3, π′ = 1/2; only q changes)
  q = (1 + ε)·q_slice, q_slice = λ√(2/3)·m = √2450894/35736;  q − q_slice = η/35736, η = ε√2450894 = η₀·√δ
  point  crossing  ε                                   η₀              δ (square-free)
  A      lower     +√(81608/27765334551) ≈ +0.0017144   +404/33736737   50234001393 = 3·1489·1907·5897
  A−     lower     −√(81608/27765334551) ≈ −0.0017144   −404/33736737   50234001393
  B      upper     +√(81608/878393760583) ≈ +0.00030480 +404/1067307121 1589220303169 = 19·41·1489·1370099
  B−     upper     −√(81608/878393760583) ≈ −0.00030480 −404/1067307121 1589220303169
  ("lower"/"upper": for that ε², the smaller/larger of the two crossings into which the parent crossing splits.)
  Exact or approximate: exact. The chain does not depend on π′ (checked at π′ = 0, 1/2, 2, −3/7 at every point).

Crossing (A and A− share ξ and r; B and B− share ξ and r)
  Covector ξ (frame components at y₀), rational:
    A, A−: ξ = (42108589454171386/227779105414740649, 186294557517875126/2050011948732665841,
                501338840/1431786279, 1371952277/1431786279)
           U = |ξ_a|²/|ξ_b|² = 13484820232/330742630449,   r = |ξ| = 344227450681/330742630449
    B, B−: ξ = (62858493798023514479/310102446766203135396, 52164981143920718591/1395461010447914109282,
                53883658045/52829177742, 591024931/26414588871)
           U = 426978392845/10460177192916,   r = |ξ| = 10887155585761/10460177192916
    (parent: U* = 400/9801, r* = 101/99)
  Eigenvalue r: rational; the 4-fold light branch meets one sheet of the other branch, so r has
    algebraic multiplicity 5 on C¹⁶ and on ker L.
  Number field: K = Q(i, √δ) with δ as above (defining polynomials x² + 1, x² − δ).

Constrained structure (on ker L(ξ)), identical at all four points
  dim ker L(ξ) = 12
  dim ker (P − r)^j on ker L, j = 1, 2, 3, 4: 4, 5, 5, 5  -> one Jordan block of size 2 (and three of size 1)
  (parent witness on the slice: 4, 6, 6 -> two blocks of size 2)
  Chain v → w with (P − r)v = 0, (P − r)w = v, v ≠ 0, Lv = Lw = 0, exact entries in K:
    data/point_A_chain.json, data/point_Am_chain.json, data/point_B_chain.json, data/point_Bm_chain.json
  Full characteristic polynomial of P(ξ): rational coefficients, (x ∓ r)⁵ (x ∓ r/2)² (x² − ω²) with ω² > 0:
    A, A−: ω² = 8432899578095448322868618474126489648820805510280575921/7784648606569962157554808331718023989682430777013573953
    B, B−: ω² = 9007066381567809668815424255239530422263007451547247402503/8314524884793793442784839403155416574079758667753171654384
    -> every characteristic speed at ξ is real. (The characteristic polynomial is the same for ±ε, as the
       step 4 dispersion relation is even in ε; the chains differ.)
  Engines: SymPy 1.14.0 (DomainMatrix over QQ<√δ + i>): scripts/step5_point_exact.py
  Cross-checked in a second engine: yes, at all four points — scripts/verify_point_fractions.py (pure Python
    fractions, own arithmetic in K, no SymPy or rs_model; rebuilds Γ, S^ρ, M, L from the paper's formulas;
    re-checks chain, det S⁰ and ranks; reproduces the parent chain on the slice; rejects tampered chains).
    Both engines are AI-written code of this project: arithmetic cross-check, not independent confirmation.

Hypotheses re-checked at these points (the paper's argument needs all four)
  det S⁰ ≠ 0: yes, exact and rational (it depends on ε² only):
    A, A−: 3487022434026159284898717563348497/150353979681895762264745351190810624
    B, B−: 43086618709359531992507710642588513/1857815457047229125650305782732230656
    each equal to the step 2 formula (9803576ε² − 763140625)(39214304ε² − 4014616321)/132101812606622109696.
  Reduction identities with q free: yes — step 2 (q symbolic), which covers every q.
  Coercivity on some ball: yes, by the paper's proof of Prop. 4.2 read with q free: its inputs are
    LL† = 3|ξ|²I and the frame E(y) (neither contains q) and m₀ = sup‖M₀‖ < ∞ (smooth coefficients, any real q).
    q changes the value of m₀ and r₀, not the argument. Checked by reading the proof, not by computation.
  Curved-symbol identification: yes, exact, at each point's ξ, with q symbolic (so for all q, in particular these).

Reproduce (repository root; Python 3, SymPy 1.14.0, mpmath)
  python3 -I scripts/points.py                                                   (constants; all points)
  python3 -I scripts/step5_point_exact.py --point A --out results/... --json data/point_A_chain.json   (~2 min; likewise A-, B, B-)
  python3 -I scripts/verify_point_fractions.py data/point_A_chain.json --out results/...              (~3 s each)
  python3 -I scripts/step5_curved_symbol.py --point A --out results/...                               (~12 s each)
  python3 -I scripts/step2_reduction_q_free.py ; python3 -I scripts/step4_dispersion_and_crossings.py
  Files, SHA-256 (run of 2026-11-22):
    1b78d24305758c2a020bfaed35b032643a485ff961d91bf3e1b811ae2399dc2b  scripts/points.py
    76cdd39e6b7ef4cb71390e6b0303e6e008b55f0f206a1c39de3e3bd2ca2310aa  scripts/rs_model.py
    8c178a5563f55a14554b677f78bf49b33f9769a28166002fe675b0ad77575bf6  scripts/step5_point_exact.py
    7aede38009faeafb08693bcf0a034d1b2e6c62302bab901905a9a047e9a932a7  scripts/verify_point_fractions.py
    5c89abd0e3ae68839ec3894aebfccc3273e48364a46209c53b9c315ed609ecf2  scripts/step5_curved_symbol.py
    4f279a6a11a638980e4e8fdfde126fe862f40e1af6ec8a83fb2c5e097dad641a  data/point_A_chain.json
    a9f1c12da5eae7ae36c39016e0a1a9d85d6d49c7286a1e017b293db59126aa8a  data/point_Am_chain.json
    67ba185bae59b6347d3959d43cba8173e969833467529df34aa4e66bbf162fa2  data/point_B_chain.json
    029596d7cb79672342cdcfc4a108da6da32717945faa5da17ba33f3df959e130  data/point_Bm_chain.json
    (run logs: results/step5_point_*_exact_, verify_point_*_fractions_, step5_point_*_curved_symbol_2026-11-22.txt;
     each log records the hashes of the scripts it ran.)
  Python 3.13.18, SymPy 1.14.0, mpmath 1.3.0, Linux x86_64. All exit codes 0.

What this does NOT establish
  - Persistence on an open set of parameters. These are four points (two values of ε², both signs, one crossing
    per point); each is a pointwise witness. No family theorem (PLAN Milestone 2 not reached).
  - Which of existence, uniqueness or the bound fails.
  - Anything about physics.
  - A re-proof of the analytic theorem at these points by anyone other than this AI-checked reading: the algebraic
    inputs are exact; coercivity is by reading the existing proof with q free.

Review status: AI-checked only. Not peer reviewed, no human expert review, not independently replicated,
not formally verified.
```

## In plain words

On the slice, the parent paper found two Jordan blocks at one crossing. Off the slice, that crossing splits in two (step 4, exact). Four points were checked exactly: both new crossings (lower and upper), each with the charge slightly above and slightly below the slice. At every one there is a Jordan block in the constraint kernel, one block where the parent had two. So the two blocks seem to separate onto the two new crossings. "Seem" because that last sentence is a pattern seen at four points and in numerics, not a proof for all small ε.

Each point was chosen so the direction ξ and the speed r are rational. The only irrational number left is one square root coming from the charge, so every check is exact arithmetic.
