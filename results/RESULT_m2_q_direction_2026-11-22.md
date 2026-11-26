# Result: Jordan chains at both off-slice crossings for every 0 < |ε| ≤ 1 (one-parameter family in q)   (date: 2026-11-22)

This is a Milestone 2 result in one direction of parameter space only: the charge q varies, and m, λ, b₁, b₂ stay at the parent witness values. For each ε it gives the same kind of pointwise statement as the parent theorem, now at every ε in an interval, not at a few points.

```
Outcome type: POSITIVE (constrained Jordan chain found) — at both crossings, for every ε with 0 < |ε| ≤ 1

Parameter family
  λ² = 1/2, m² = 823/142944, b₁ = 1, b₂ = 1/3, any π′;  q = (1 + ε)·q_slice, q_slice = λ√(2/3)·m,  0 < |ε| ≤ 1.
  (ε = 0 is the parent witness.)

Statement (algebraic part)
  For every such ε, the off-slice characteristic crossing splits into two real crossings, at covectors
  ξ±(ε) = (s±(ε), 0, 1, 0) with s±(ε)² = u±(ε) > 0 the two roots of the crossing polynomial C(u, ε²). The same
  holds at every covector obtained from these by rotations within each sphere factor and by scaling. At each of
  them, r = |ξ| is an eigenvalue of P(ξ) restricted to ker L(ξ) with algebraic multiplicity ≥ 5 and geometric
  multiplicity ≤ 4. So P|ker L is not diagonalisable at r: there is a nontrivial Jordan chain v → w in ker L(ξ),
  with v ≠ 0, (P − r)v = 0, (P − r)w = v. It is the same chain for every π′.

Proof structure (each piece exact; scripts and logs in this repository)
  (R) Rotation equivariance: P(Rᵀξ)Ŝ = ŜP(ξ), L(Rᵀξ)Ŝ = S L(ξ), Ŝ = Rᵀ ⊗ S, S = cos(θ/2) + sin(θ/2)ΓⁱΓʲ, for rotations
      in the (1,2) and (3,4) planes; symbolic ξ, q, θ.  scripts/m2_rotation_invariance.py
  (M) Algebraic multiplicity ≥ 5: along the ray, charpoly(P|ker L) ∝ (x² − s² − 1)⁴·Q(x², s², ε²) exactly, with ε
      symbolic; at a crossing, Q(r², s², ε²) = 0 adds one root r to the four light-branch roots.
      scripts/step4_dispersion_and_crossings.py
  (C) Crossings real, distinct, positive for 0 < ε² ≤ 1: the discriminant of C in u is
      95227628033791228445188816896·ε²(91809 − 823ε²) > 0, the u² and u⁰ coefficients are > 0 and the u¹
      coefficient < 0 on [0, 1] (exact Sturm counts).  scripts/m2_family_q_direction.py
  (G) Geometric multiplicity ≤ 4. Parametrise s = 2τ/(1 − τ²), so |ξ| = (1 + τ²)/(1 − τ²); write η = ε√2450894.
      The polynomial matrix Dn·(P − |ξ|)K(τ, η) (K a basis of ker L; Dn = (τ² − 1)²(4η² − 27625²)(16η² − 63361²),
      nonzero on the region used) has a fixed 8×8 minor d. With τ = 1/10 + a (τ = 1/10 ↔ s* = 20/99) and η = b,
      d(a, b) has degree (32, 32), 561 terms and d(0, 0) > 0. On the box |a| ≤ 1/20, |b| ≤ 782767/500 (≥ √2450894):
      Σ_{(i,j)≠0} |c_ij| hₐ^i h_b^j < |d(0, 0)| (exact rational inequality, ratio 0.899). So d ≠ 0 on the box and
      rank (P − r)|ker L ≥ 8 there.  scripts/m2_family_q_direction.py, data/m2_minor.json
  (B) Both crossings stay in the box for |ε| ≤ 1: the crossing polynomial has no real root on the box sides
      τ = 1/10 ± 1/20 for |η| ≤ h_b (exact Sturm count), and the two crossing branches are continuous in ε and
      start at τ = 1/10. (Smaller boxes certify the same for |ε| ≤ 1/1000, 1/100, 1/10, 1/2.)
  Together: geometric ≤ 4 < 5 ≤ algebraic at every crossing for 0 < |ε| ≤ 1.

Hypotheses of the paper's argument, for every such ε (each ε gives its own pointwise theorem)
  det S⁰ ≠ 0: yes — det S⁰(ε) = (9803576ε² − 763140625)(39214304ε² − 4014616321)/132101812606622109696 ≠ 0 for
     |ε| < 8.82 (step 2).
  Reduction identities with q free: yes — step 2 (q symbolic).
  Curved-symbol identification: yes, exactly, for EVERY chart covector k and EVERY q (k and η symbolic):
     scripts/m2_curved_symbol_general.py.
  Coercivity on some ball: by the paper's proof of Prop. 4.2 read with q free (inputs: LL† = 3|ξ|²I and the frame,
     neither of which contains q, and m₀ = sup‖M₀‖ < ∞ for any real q). Proof reading, not computation. The radius
     r₀ may depend on ε; no uniformity in ε is claimed.

Engines and cross-checks
  SymPy 1.14.0 for (R), (M), (C), (G), (B) and the curved symbol.
  Second engine (pure Python fractions, no SymPy): scripts/verify_m2_fractions.py
    - rebuilds S^ρ, M, B(ξ) from the paper's formulas and evaluates the 8×8 minor at a 33 × 33 grid of rational
      (τ, η). It equals the recorded d(a, b) at all 1089 points. Both have degree ≤ 32 in each variable, so they are
      the same polynomial;
    - re-checks the box inequalities and the box-side root counts with its own Sturm sequences;
    - negative controls: a tampered coefficient and a too-narrow box both FAIL.
  (M) the ray dispersion relation: since cross-checked in the fractions engine (verify_m2_dispersion_fractions.py:
    the full 16×16 characteristic polynomial along the ray, proved on a degree-sufficient 17 × 65 grid; added later on
    2026-11-22, see RESULT_m2_extended).
  Not cross-checked in a second engine: (R) rotation equivariance and the general curved-symbol identity (SymPy only). The four exact points A, A−, B, B− (RESULT_offslice_points) are
    special cases, checked in both engines.
  Both engines are AI-written code of this project: arithmetic cross-check, not independent confirmation.

Reproduce (repository root; Python 3, SymPy 1.14.0)
  python3 -I scripts/m2_rotation_invariance.py --out results/...            (~10 s)
  python3 -I scripts/m2_family_q_direction.py --out results/... --json data/m2_minor.json   (~30 s)
  python3 -I scripts/verify_m2_fractions.py data/m2_minor.json --out results/...           (~2 min)
  python3 -I scripts/m2_curved_symbol_general.py --out results/...          (~8 s)
  python3 -I scripts/step4_dispersion_and_crossings.py ; python3 -I scripts/step2_reduction_q_free.py
  Files, SHA-256 (run of 2026-11-22):
    6533e699c5dbab44cb59d15befdb2ec0a031a0de3276075c0d4a1582f7a8d119  scripts/m2_rotation_invariance.py
    c38c85bb424712403b7250cf651ec3e077be23560f86540f2652d77f7dcbfa11  scripts/m2_family_q_direction.py
    dff71899ff256efcced35ce7b8139fa9905345bb899fc20de0d2848c2a04f0d9  scripts/verify_m2_fractions.py
    d9b8d320cba54d2978609b2c2042073160cb73be2697445e2f7a2fa3b0d902de  scripts/m2_curved_symbol_general.py
    faa705683310ca2f9d86885e8072b6ddab332f14ebb79b07bae5a19edfdfbd69  data/m2_minor.json
    77f43145d64c0e4fe3894fd5ec45b60abd4f15bd647ef85bf36fab80bc462d49  results/m2_rotation_invariance_2026-11-22.txt
    420ece35a99697d48712eb7cb127190cfd18588a3df38e96f75641201a6cf43c  results/m2_family_q_direction_2026-11-22.txt
    a24f1308cb638ff6c011bb8f8f43da8cd7de5c66878122a0413f0abb8545d7eb  results/verify_m2_fractions_2026-11-22.txt
    3462b290951e9f647a10b4b812670dc2baf06c9a6965746e3c708d665be6d963  results/m2_curved_symbol_general_2026-11-22.txt
  Python 3.13.18, SymPy 1.14.0, Linux x86_64. All exit codes 0.

What this does NOT establish
  - Anything off this one-parameter curve: m, λ, b₁, b₂ were not varied. Not an open set in parameter space.
  - Uniformity in ε of the radius or constants: each ε is a separate pointwise theorem.
  - |ε| > 1 (not attempted; det S⁰ ≠ 0 holds up to 8.82, the crossings stay real up to 10.56).
  - Which of existence, uniqueness or the bound fails; anything about physics.
  - A re-proof of the analytic theorem by anyone other than this AI-checked reading (coercivity by proof reading).

Review status: AI-checked only. Not peer reviewed, no human expert review, not independently replicated,
not formally verified.
```

## In plain words

Change the charge by any amount up to ±100% of its slice value (ε between −1 and 1), keeping everything else as in the paper. The paper's crossing then always splits into two. At both of them, the symbol restricted to the constraint kernel has a Jordan block. That is the ingredient the paper's argument needs, and the other ingredients were re-checked for every such charge. The four exact points of the earlier result are special cases.

The proof does not compute anything at each ε separately. It shows that one fixed determinant cannot vanish anywhere in a box that provably contains all the crossings. That is a finite, exact inequality between rational numbers, and two independent arithmetic engines check it.
