# Partial identification of classifier dependence from labeled accuracies + unlabeled data: sharp bounds and robust weights for weighted votes

Conventions: Z_i = 1{classifier i correct}, π = law of Z on {0,1}^K, p_i = P(Z_i=1), f = prediction vector, Y ∈ {0,1}, z̄ = 1−z (antipodal pattern). φ(v) = P(f=v) is the law of the predictions, which unlabeled data give. α(z)=P(Y=1,Z=z), β(z)=P(Y=0,Z=z), η = P(Y=1). A rule g:{0,1}^K→{0,1} is *self-dual* if g(v̄)=1−g(v). Every weighted vote with positive weights, no bias and no ties is self-dual. G = 2g−1, s = 2z−1.

All code is in `S = code/unlabeled_thread/`. The library is `core.py`. Experiments are `q1_identification.py`, `q1b_walsh_check.py`, `q2_robust_weights.py`, `q2b_regret_and_anticorr.py`, `q2c_phase_symmetric.py`, `q3_ml_bounds.py`, `q3b_closed_forms.py`, `q4_empirical.py`, `q4b_identification_only.py`, `q4c_p_only_conf.py`, `count_selfdual_monotone.py` and `make_figures.py`. Raw outputs are in `S/results/*.json|npy|txt`. Figures were copied to this folder as `unlabeled_fig1_sharp_vs_pacbayes.png` … `unlabeled_fig5_empirical.png`. All LPs use scipy 1.18.1 `linprog(method="highs")`, numpy 2.5.3 and scikit-learn 1.9.1. Labels used below: **[PROVED]** means a short proof is written here; **[NUM]** means numerically verified, with the instance count and maximum error given; **[SPEC]** means speculative.

---

## Q1. What unlabeled data + marginals identify about π, and the sharp interval for a weighted vote

### Takeaway
With the class prior unknown, the law of f plus the marginals p identify π only up to the polytope **E(m,p) = {π ≥ 0 : π(z)+π(z̄) = m(z), Σ_z z_i π(z) = p_i}**, where m(z) = φ(z)+φ(z̄) is the "even part". A known prior adds two total-variation constraints, and these rarely shrink the interval much. The candidate claim is confirmed with one important refinement. When Z ⟂ Y and η ≠ ½, π is fully identified, so "only the even part" holds in case (b) only for balanced classes. The third-order success t = P(Z1=Z2=Z3=1) is not identified in case (a): the interval had positive width in 500/500 random instances. For K=3 the full unlabeled law carries exactly the same information as the three pairwise agreements, so the identified set is the known segment t ∈ [L,U]. For K=5 the full unlabeled law cuts the width of the majority-accuracy interval roughly in half compared with pairwise agreement alone.

### Cited Findings
**Basic identities [PROVED].**
- If Y=1 then f=Z, and if Y=0 then f=Z̄. Hence φ(v) = α(v)+β(v̄), and φ(v)+φ(v̄) = π(v)+π(v̄) =: m(v). The pairwise agreement P(f_i=f_j) = P(Z_i=Z_j) is an even function of z, so it is a linear functional of m. This gives r_ij = (P(f_i=f_j) − 1 + p_i + p_j)/2. — [core.py](S/core.py)
- For self-dual g, 1{g(f)=Y} = g(Z) holds pointwise. Hence acc(g) = Σ_z π(z)g(z) = ½ + ½ Σ_z [π(z)−π(z̄)] g(z). The proof uses Σ_z ½(π(z)+π(z̄))g(z) = ½, which follows from g(z̄)=1−g(z). The candidate formula is confirmed. In Walsh form, acc(g) = ½ + ½ Σ_{|S| odd} Ĝ(S) E[χ_S(s)], because Ĝ vanishes on even degrees. Numerical check of the identity: max error 2.2e−16 over 200 instances each for K=3 and K=5. — [q1b_walsh_check.py](S/q1b_walsh_check.py), `results/q1b_walsh.json`

**Theorem 1 (prior unknown, case (a)) [PROVED + NUM].** The set of π compatible with (φ, p) is exactly E(m,p).
- *Proof.* Fix an antipodal pair {z,z̄} and π in E(m,p). The free unknowns are α(z) ∈ [0,π(z)] and α(z̄) ∈ [0,π(z̄)], with β = π − α. The two φ-equations reduce to one condition, α(z) − α(z̄) = φ(z) − π(z̄). The left side ranges over [−π(z̄), π(z)], and the condition is solvable iff 0 ≤ φ(z) ≤ m(z), which always holds. ∎
- *Numerics.* A direct LP over (α, β) ∈ R^{2·2^K} with constraints φ and p was compared with the LP over E(m,p), on random linear objectives. Max |difference| was 1.1e−15 (K=3), 1.1e−15 (K=4) and 2.0e−15 (K=5), over 120/120/60 random joint laws with class-dependent error patterns × 4 objectives. — [q1_identification.py](S/q1_identification.py), `results/q1_results.json`
- *Equivalent "pair-splitting LP".* Use variables x_v = π(v) ∈ [0, m(v)] for v in one half of the cube (2^{K−1} variables) and set π(v̄) = m(v) − x_v. The K marginal constraints read Σ_v x_v(2v_i−1) = p_i − Σ_v m(v)(1−v_i). The accuracy of a self-dual g is const + Σ_v x_v·G(v). This is a box-constrained LP with K equalities, so a basic optimal solution has at most K pairs strictly inside their box. **Structure [PROVED]:** the worst-case adversary "orients" every antipodal pair, putting all mass on z or on z̄, except at most K pairs. — [core.py `set_even`](S/core.py)
- *Dimension [NUM].* dim E(m,p) = 2^{K−1} − K, which is 1 for K=3 and 11 for K=5. dim P(p,r) = 2^K − 1 − K − C(K,2), which is 1 for K=3 and 16 for K=5 (matrix-rank check). — `results/q1b_walsh.json`

**Theorem 2 (which functionals are point-identified) [PROVED].** If π has full support, the direction space of E(m,p) is {δ odd, δ ⊥ s_i ∀i}. This is the span of odd Walsh characters of degree ≥ 3. So a linear functional c·π is identified iff c has no Walsh weight on odd degrees ≥ 3. Corollaries:
- Among self-dual rules, only the K dictators have point-identified accuracy. Majority accuracy is never point-identified when π has full support.
- Under pairwise information P(p,r), the direction space is all characters of degree ≥ 3. The full unlabeled law therefore adds exactly the even-degree moments of degree ≥ 4, which exist only for K ≥ 4. For K=3 the two sets coincide: E(m,p) = P(p,r), with max LP difference 6.7e−16 over 3000 instances. — `results/q1_results.json` (key C_K3_closed_form)

**Theorem 3 (prior known, case (a)) [PROVED + NUM].** The identified set is E(m,p) ∩ {TV(π, φ) ≤ 1−η} ∩ {TV(π, φ∘flip) ≤ η}.
- *Proof.* On each pair, the achievable α-mass α(z)+α(z̄) ranges over the interval [|φ(z)−π(z̄)|, m(z) − |π(z)−φ(z)|]. Summing over pairs gives η ∈ [TV(π, φ̃), 1 − TV(π, φ)]. ∎
- *Numerics.* The LP with these TV constraints matched the direct (α,β,η) LP to max 1.7e−15 (K=3,4,5).
- *Value of knowing η.* The mean ratio of interval widths, known η over unknown η, for random linear functionals was 0.970 (K=3), 0.943 (K=4) and 0.959 (K=5), with medians 1.000/1.000/0.999. — `results/q1_results.json` (AB_K*)

**Third-order t (K=3, case (a)) [NUM].** The interval for t had positive width in 500/500 random instances whether η was unknown or known. Mean width was 0.0784 with η unknown and 0.0751 with η known. Knowing η shrank the interval in only 20.2% of instances. — `results/q1_results.json` (D_t_identification)

**Case (b), Z ⟂ Y [PROVED + NUM].**
- Here φ = ηπ + (1−η)π̃, so φ − φ̃ = (2η−1)(π − π̃). Also 2P(f_i=1) − 1 = (2η−1)(2p_i−1), so η is identified from any p_i ≠ ½.
- If η ≠ ½, π is fully identified as π = ½[m + (φ−φ̃)/(2η−1)]. The recovery error was max 1.5e−14 over 400 instances (K=3, 5). This matches the known "independent errors ⇒ accuracies identifiable from agreement rates" phenomenon. — `results/q1_results.json`; [Platanios, Blum & Mitchell, UAI 2014](https://proceedings.mlr.press/r12/university14g.html)
- If η = ½ (balanced classes with symmetric class-conditional errors), φ = m/2 and the identified set is again E(m,p). The candidate "even part only" statement is exactly right in this sub-case.

**K=3 closed form in observables [PROVED + NUM].**
- Let D = P(f not unanimous) and d_i = P(classifier i is the lone dissenter) = (D_ij + D_ik − D_jk)/2, where D_ij = P(f_i ≠ f_j). All of these come from unlabeled data.
- Majority accuracy A is affine in the single free parameter. The sharp interval is **A ∈ [max{max_i(p_i − d_i), Σp − 2 + D}, min{min_i(p_i + d_i), Σp − D}]**.
- *Proof.* Write x_0 = π(111) and x_i = π(only i wrong). Then x_i = (A + d_i − p_i)/2 ∈ [0,d_i] and x_0 = (Σp − D − A)/2 ∈ [0, 1−D]. ∎
- This matched the LP to max 1.0e−15 over 3000 instances and is the observable form of the known [R−2U, R−2L]. — `results/q1_results.json`

**Dual "unlabeled-margin" certificate [PROVED by LP duality; NUM].**
- For self-dual g, the sharp lower bound is acc_min(g) = max_{λ∈R^K} ½[1 + λ·(2p−1) − E_{f∼φ}|1 − G(f)·λ·s(f)|]. This matched the primal to max 1.8e−15 over 250 instances.
- For a weighted vote with weights w, restricting to λ = c·w gives the bound acc ≥ ½[1 + c w·(2p−1) − E_f|1 − c|w·s(f)||] for every c ≥ 0. It needs only p and the unlabeled margin distribution |w·s(f)|. It is valid but not sharp: the gap to the sharp bound had mean 0.061, median 0.050 and max 0.228, and it was exact in 0/250 instances.
- This is the population, equality-constraint analogue of Lemma 7 of Balsubramani & Freund, min_z (1/n)z·g = max_{σ≥0}[b·σ − (1/n)‖Fᵀσ − g‖₁]. — `results/q1_results.json` (F_duality); [Balsubramani & Freund, COLT 2015](https://arxiv.org/pdf/1503.01811)

**K=5 value of the full unlabeled law [NUM].** Over 300 random class-dependent joint laws, the mean width of the 5-voter majority interval was:

| Information | Mean width |
|---|---|
| p + pairwise agreement | 0.2867 |
| p + full unlabeled law | 0.1606 (mean ratio 0.538 vs pairwise) |
| full unlabeled law + known η | 0.1503 |

The mean gain in the lower bound from full vs pairwise information was 0.0636. The truth was inside every interval. — `results/q1_results.json` (G_K5…), Fig. `unlabeled_fig4_K5_intervals.png`

### Inferences
- The "even part / odd part" picture is the correct structural lens. Unlabeled data pin the even Walsh moments, labels (p) pin the degree-1 odd moments, and accuracy depends only on odd moments. Partial identification enters only through positivity, |π(z) − π(z̄)| ≤ m(z), which links the odd part to the identified even part. That is why the bound is a box LP ("pair splitting") and not vacuous.
- Knowing η is almost worthless in case (a) (median width ratio ≈ 1). This suggests reporting results for the prior-unknown set, which is also simpler.
- Theorem 2 gives a crisp negative message for a paper: with labels used only through p and any amount of unlabeled data, no non-dictator vote has an identified accuracy. Every ensemble claim is inherently an interval claim.

### Gaps
- For general K, I did not find a closed form for the pair-splitting LP. I only have the "≤ K fractional pairs" structure and the dual formula. Whether weighted votes admit a sorting/greedy solution, as in fractional knapsack, is open.
- Identification when labels give more than p, for example class-conditional accuracies P(Z_i=1|Y=y), was not analysed. That case breaks the self-dual symmetry.

---

## Q2. Robust (max–min) weights from (p, pairwise agreement), and comparison with log-odds weights

### Takeaway
For K=3 the max–min rule has an exact closed form: **use majority iff Σ_{i<j} T_ij < min_i e_i, otherwise the best single classifier**. Here T_ij = P(both i and j err) = (e_i + e_j − D_ij)/2 comes from labeled errors plus unlabeled disagreement. Mixing or randomising rules never helps for K=3. For K=5 the max–min rule over the pairwise set, on positively dependent random instances, was always either a dictator or a 3-voter sub-majority. It differs from log-odds weights in about 94% of instances. It raises worst-case accuracy by about 7 points on average and costs about 1 point of average true accuracy on these generators. It wins outright on negatively dependent ("disjoint error regions") instances. The pairwise identified set is wide, so max–min is very conservative. Minimax regret is a less conservative alternative and sometimes beats log-odds (bloc family).

### Cited Findings
**Rule enumeration [NUM].** Brute force shows every self-dual monotone Boolean function on K ≤ 5 variables is a positive-weight threshold function. The counts are 4 (K=3), 12 (K=4) and 81 (K=5); the K=5 rules fall into 7 weight types. The random-weight plus integer-weight enumeration found exactly these. — [count_selfdual_monotone.py](S/count_selfdual_monotone.py)

**Theorem 4 (K=3 max–min rule) [PROVED + NUM].**
- Dictator i has identified accuracy p_i. Majority's worst case over P(p,r) is max{max_i(p_i−d_i), Σp−2+D} = max{…, 1 − Σ_{i<j}T_ij}, using the identity D = Σe − ΣT, which holds to max error 6.7e−16. Majority errs iff some pair errs jointly, so 1 − ΣT is the tight Bonferroni/union bound.
- Hence max–min picks majority iff ΣT_ij < e_min, equivalently D > e_(2) + e_(3), with the errors sorted ascending. The max–min value is max{p_max, 1 − ΣT}.
- The rule matched the LP in 834/834 (q2) and 2000/2000 (q3b) instances.
- *Symmetric case:* majority is optimal iff ρ < ρ*(p) = ((5p−2)/3 − p²)/(p(1−p)). Example values: ρ*(0.6) = −0.111, ρ*(0.7) = 0.048, ρ*(0.8) = 0.167, ρ*(0.9) = 0.259, ρ*(0.95) = 0.298, with ρ* → 1/3 as p → 1. At ρ=0 the rule picks majority iff p > 2/3. The threshold agreed with the LP on the whole (p,ρ) grid (agreement 1.0000). — `results/q2_results.json`, `results/q3b_results.json`, `results/q2c_results.json`, Fig. `unlabeled_fig2_phase_symmetric.png`
- *No gain from randomisation for K=3 [PROVED].* Dictator accuracies are constant on the 1-D identified set and majority is affine in t, so all worst cases sit at the same endpoint. Numerically, mixtures of rules, randomised self-dual g and the Balsubramani–Freund-type game over all randomised predictors gave gain 0.0000 (max) over 834 instances.

**Random-instance comparison (q2; families: independent / latent-mixture / "bloc" of correlated voters; instances with some p_i ≤ ½ discarded) [NUM].** "worst" is the worst case over P(p,r); "true" is the accuracy under the generating law. — `results/q2_results.json`, [q2_robust_weights.py](S/q2_robust_weights.py)

| K | family (n) | max–min ≠ log-odds | true: max–min / log-odds / majority / best single | worst: max–min / log-odds / majority / best single |
|---|---|---|---|---|
| 3 | independent (300) | 0.320 | 0.8318 / 0.8376 / 0.8249 / 0.8183 | 0.8262 / 0.8162 / 0.7871 / 0.8183 |
| 3 | mixture (234) | 0.376 | 0.8309 / 0.8357 / 0.8177 / 0.8202 | 0.8256 / 0.8096 / 0.7703 / 0.8202 |
| 3 | bloc (300) | 0.773 | 0.7628 / 0.7738 / 0.7629 / 0.7558 | 0.7585 / 0.7051 / 0.6852 / 0.7558 |
| 5 | independent (200) | 0.885 | 0.8836 / 0.8959 / 0.8756 / 0.8464 | 0.8723 / 0.8168 / 0.7591 / 0.8464 |
| 5 | mixture (143) | 0.944 | 0.8842 / 0.8936 / 0.8726 / 0.8547 | 0.8739 / 0.8136 / 0.7639 / 0.8547 |
| 5 | bloc (200) | 0.985 | 0.8034 / 0.8140 / 0.7972 / 0.7831 | 0.7930 / 0.6921 / 0.6624 / 0.7831 |

- *Direction of disagreement.* For K=3, all 416 disagreements in these families had log-odds choosing majority and max–min choosing a dictator, with 0 the other way.
- *When the two differ.* The max–min rule had higher true accuracy in only 15.9% (K=3) and 14.2% (K=5) of disagreements. Its worst-case gain averaged +0.055 (K=3) and +0.078 (K=5).
- *Rule types chosen, K=5.* Max–min chose only (1,1,1,0,0) in 315 instances and (1,0,0,0,0) in 228. Log-odds spread over all types; the most common was simple majority, 165/543.
- *Randomisation, K=5.* A mixture of rules gained mean 0.0009 (max 0.0317) over the best deterministic rule on P(p,r). The BF game on the full unlabeled set gained mean 0.0050 (max 0.0478). Using the full unlabeled law, the max–min rule differed from the pairwise one in 21.9% of instances and raised mean true accuracy from 0.8542 to 0.8564.
- *Sub-majority worst case [NUM].* For K=5, the worst case of a 3-voter sub-majority over the 5-voter P(p,r) equals the K=3 closed form on that triple (200/200, max difference 4.4e−16). The worst case depends only on the triple's own (p, r). — `results/q3b_results.json`

**Negatively dependent family and minimax regret (q2b; independent RNG, so numbers differ slightly from q2) [NUM].** — `results/q2b_results.json`, `results/q2b_log.txt`
- *"anti" family:* errors occupy arcs of a circle, so error sets are nearly disjoint.
  - K=3 (300): max–min 0.8856 = oracle, log-odds 0.8639, majority 0.8553, best single 0.8252. Max–min > log-odds in 37.3% of instances and < in 0%. In 13.7% max–min chose majority while log-odds chose a dictator.
  - K=5 (150): max–min 0.9245 = oracle, log-odds 0.8791. Max–min > log-odds in 73.3% and < in 0%.
  - *Caveat:* this family is degenerate. The majority interval had **zero width in 200/200 instances for both K=3 and K=5**, because the true law sits at a vertex of the identified set (for example, no sample where all err). — `results/q4c_summary.json`
- *Minimax-regret mixtures (true accuracy):*
  - K=5 independent: regret-mixture 0.8767, regret-det 0.8786, log-odds 0.8861, max–min 0.8743.
  - K=5 mixture: regret-mixture 0.8850, regret-det 0.8874, log-odds 0.8919, max–min 0.8807.
  - K=5 bloc: regret-mixture 0.8103, regret-det 0.8122, log-odds 0.8070, max–min 0.7981.
  - The minimax regret value over P(p,r) is large: 0.0837, 0.0751 and 0.1049 for K=5 in these three families, and 0.015–0.024 for K=3.

**Symmetric phase diagram (p ∈ [0.505, 0.995], ρ ∈ [−0.5, 1]; one representative per weight type) [NUM].** — `results/q2c_results.json`
- K=3: the dictator is optimal on 72.6% of feasible cells and majority on 27.4%.
- K=5: dictator 77.3%, 3-sub-majority 17.9%, full 5-majority 4.8%. The 5-majority wins only in a band of negative ρ (for example ρ ∈ [−0.25, −0.09] at p=0.6). At ρ = 0 the 3-sub-majority is optimal for p ≳ 0.70 and the dictator below.
- Log-odds weights are equal in this case, so log-odds equals majority everywhere.

### Inferences
- Robust weights given only (p, pairwise) are *sparse*. They concentrate on the smallest sub-committee whose worst case can be certified from the information available. A plausible explanation [SPEC] is Theorem 2: a rule's accuracy splits into an identified degree-1 part and a non-identified odd part of degree ≥ 3. Small committees (one voter or three) have fewer unidentified high-order terms, and for triples the worst case is exactly computable from pairwise data.
- Log-odds weights are optimal only under independence, and their worst case over the pairwise set is poor: about 0.08–0.11 below max–min for K=5 positive-dependence families. Max–min, however, gives up about 1 point of average accuracy under benign dependence. This is the classic robustness/efficiency trade-off. Minimax regret, or max–min over a set narrowed by extra information (full unlabeled law, priors on higher-order dependence), is the natural middle ground and a candidate paper contribution.
- The K=3 rule "ΣT < e_min" is simple and interpretable. It needs labeled errors plus unlabeled disagreement, and it reads naturally as "majority is safe iff joint failures are rarer than the best voter's failures".

### Gaps
- I have no closed form for the K=5 max–min rule in general. The "dictator or triple" observation holds on 543/543 positively dependent random instances but not in the anti family or the negative-ρ symmetric band.
- I did not study finite-sample-robust max–min with confidence sets in synthetic settings; Q4 does this only on real data.

---

## Q3. Sharp vs existing PAC-Bayes majority-vote bounds

### Takeaway
Let μ = E_ρ[L(h)] and s = E_{ρ²}[L(h,h')] (tandem loss, diagonal included). As K → ∞ the sharp bound is **L(MV) ≤ min{C-bound, 3μ − 2s}**, and 3μ − 2s = E_ρ[L] + E_{ρ²}[D(h,h')] (Gibbs risk + expected disagreement).
- The C-bound is sharp exactly when s ≤ μ/2, equivalently d_Q ≥ R(G_Q), which is consistent with Germain et al. Prop. 13.
- In the complementary region the sharp value is the *linear* bound μ + E[D], strictly below 2μ, 4s and C.
- For finite K the sharp bound is lower again and has an explicit Boros–Prékopa-type closed form. Over a grid of (μ,s), the best existing bound min{2μ, 4s, C} exceeds the sharp bound by a mean of 0.222 (K=3), 0.168 (K=5) and 0.095 (K→∞).
- With individual (p, r) instead of averages, the bounds tighten much further. Mean bound for K=5 uniform votes, true error 0.132: C-bound 0.423 → sharp|(μ,s) 0.247 → sharp|(p,r) 0.244 → sharp|(p, full unlabeled) 0.184.

### Cited Findings
**Existing bounds (oracle forms).**
- First-order: L(MV_ρ) ≤ 2E_ρ[L(h)] (Wu et al. 2021, Thm 1). Tandem: L(MV_ρ) ≤ 4E_{ρ²}[L(h,h')] (Thm 3, from Masegosa et al. 2020). — [Wu et al., NeurIPS 2021, arXiv:2106.13624](https://arxiv.org/pdf/2106.13624.pdf)
- C-bound: L(MV) ≤ (E_{ρ²}[L(h,h')] − E_ρ[L]²)/(¼ + E_{ρ²}[L(h,h')] − E_ρ[L]), valid if E_ρ[L] ≤ ½ (Thm 5). Chebyshev–Cantelli parametric bound: L(MV) ≤ (E_{ρ²}[L(h,h')] − 2νE_ρ[L] + ν²)/(½−ν)² for ν < ½ (Thm 7). Its minimiser ν* = μ − (s − μ²)/(½ − μ) recovers the C-bound and ν = 0 gives the tandem bound. — [Wu et al. 2021](https://arxiv.org/pdf/2106.13624.pdf)
- I checked the Thm 7 claims numerically: min_ν Thm 7 = C-bound to max 1.6e−15 over 400 points. A naive grid restricted to ν ≥ −5 misses ν*, which can be very negative. — `results/q3_results.json`
- Germain et al. state the C-bound as R(B_Q) ≤ 1 − (1−2R(G_Q))²/(1−2d_Q) = 1 − μ₁²/μ₂ (Thm 11), with ties counted as errors (R(B_Q) = P(M_Q ≤ 0)). Prop. 13 says the C-bound is attained by some margin variable with the same two moments **iff 0 < μ₂ ≤ μ₁** (equivalently R(G_Q) ≤ d_Q). Prop. 14 says this is equivalent to C_Q ≤ 2R(G_Q). The disagreement d_Q is label-free. — [Germain et al., JMLR 2015, arXiv:1503.08329](https://ar5iv.arxiv.org/html/1503.08329)
- The tandem loss counts an error only if both h and h′ err. ρ² is the product distribution, so the diagonal h=h′ is included. For binary tasks E_{ρ²}[L(h,h')] = E_ρ[L(h)] − ½E_{ρ²}[D(h,h')] (attributed to Lacasse et al. 2007). — [Masegosa et al., NeurIPS 2020, arXiv:2007.13532](https://ar5iv.labs.arxiv.org/html/2007.13532)

**Theorem 5 (sharp bound from (μ, s), continuum / K→∞, ties counted as errors) [PROVED + NUM].**
- The bound is sup{P(W ≥ ½) : W ∈ [0,1], EW = μ, EW² = s} = (s−μ²)/(s−μ+¼) if s ≤ μ/2, and 3μ − 2s if s ≥ μ/2. Equivalently it is min{C-bound, 3μ − 2s}.
- *Proof of the new branch.* The quadratic q(w) = 3w − 2w² satisfies q ≥ 0 on [0,½] and q ≥ 1 on [½,1], since q−1 = −(2w−1)(w−1). So L(MV) ≤ 3μ − 2s always. Equality holds for W supported on {0, ½, 1}, which is feasible iff s ≥ μ/2. The other branch is Cantelli, sharp by Germain Prop. 13. ∎
- *Numerics.* Matched a 2001-point grid LP to max 2.4e−6 (grid discretisation) over 400 points.
- Substituting s = μ − ½E[D] gives **L(MV) ≤ E_ρ[L(h)] + E_{ρ²}[D(h,h')]**. This beats 2E_ρ[L] exactly when E[D] < E[L]. It is linear in the moments, so it is directly estimable, with E[D] from unlabeled data. — [q3_ml_bounds.py](S/q3_ml_bounds.py), [q3b_closed_forms.py](S/q3b_closed_forms.py)

**Theorem 6 (finite odd K, uniform weights) [PROVED validity via dual quadratics; sharpness NUM].**
- Let ν = number of errors and r = (K+1)/2. Then L(MV) = P(ν ≥ r) ≤ min{1, E[ν(r+K−ν)]/(rK), min_{j=0..r−2} E[(ν−j)(ν−j−1)]/((r−j)(r−j−1))}, with E ν = Kμ and E ν² = K²s.
- Matched the LP over the distribution of ν to max 1.3e−7 (K=3), 4.3e−10 (K=5) and ≤ 1.7e−15 (K = 7, 11, 21, 101), over 1500 random moment pairs per K. — `results/q3b_results.json`

**Grid comparison [NUM].** The grid is μ ∈ {0.01,…,0.49} (49 values) × s = μ² + f(μ−μ²) for 39 interior values of f, keeping K-feasible points. "Best existing" is min{2μ, 4s, C, 1}. — `results/q3_results.json`, Fig. `unlabeled_fig1_sharp_vs_pacbayes.png`

| K | points | mean sharp | mean best existing | mean gap (max) | mean ratio sharp/best (min) | fraction where best existing is sharp |
|---|---|---|---|---|---|---|
| 3 | 1643 | 0.2877 | 0.5100 | 0.2222 (0.5504) | 0.534 (0.001) | 0.000 |
| 5 | 1819 | 0.3134 | 0.4817 | 0.1683 (0.5407) | 0.619 (0.006) | 0.000 |
| 11 | 1900 | 0.3400 | 0.4655 | 0.1255 (0.4785) | 0.725 (0.033) | 0.000 |
| 21 | 1909 | 0.3529 | 0.4633 | 0.1104 (0.4781) | 0.776 (0.232) | 0.000 |
| 101 | 1911 | 0.3645 | 0.4629 | 0.0983 (0.4776) | 0.819 (0.513) | 0.000 |
| ∞ | 1911 | 0.3676 | 0.4629 | 0.0953 (0.4775) | 0.829 (0.513) | 0.304 (= share of grid with s ≤ μ/2, 0.3035) |

Illustration at μ = 0.2:
- s = 0.08: sharp is 0.0600 (K=3), 0.1667 (K=5), 0.2493 (K=11) and 0.3077 (K=∞). The C-bound is 0.3077, 2μ is 0.400 and 4s is 0.320.
- s = 0.14: sharp is 0.2900 (K=3), 0.3000 (K=5) and 0.3200 (K=∞). The C-bound is 0.5263 and 2μ is 0.400.

**Weighted majority vote [PROVED reduction + NUM].**
- With weights ρ, W = Σ_i ρ_i(1−Z_i) takes values in V_ρ = {Σ_{i∈S} ρ_i}. The 2^K-variable LP (max P(W > ½) s.t. E_ρ moments) collapses to an LP over V_ρ: any law on V_ρ is realisable on {0,1}^K, and the objective and constraints depend only on W.
- Closed form, with w* = min{v ∈ V: v > ½}: min{1, ((w*+1)μ − s)/w*, min over consecutive a < b in V ∩ [0,½] of (s − (a+b)μ + ab)/((w*−a)(w*−b))}.
- Matched the 2^K LP to max 1.0e−15 (K=3) and 1.7e−15 (K=5) over 300 random weightings each. The continuum bound exceeds the finite-V sharp bound by a mean of 0.052 (K=3) and 0.037 (K=5).
- For K=3 with *individual* pairwise info, the sharp error bound is **min{Σ_{i<j}T_ij, min_i(e_i + d_i)}**, matching the LP to max 3.3e−16 over 2000 instances. — `results/q3b_results.json`, `results/q3_results.json`

**Information hierarchy (K=5, 287 random latent-mixture instances, oracle quantities, bounds capped at 1) [NUM].** The chain true ≤ sharp|full ≤ sharp|pairwise ≤ sharp|(μ,s) ≤ C held in every instance, and sharp|(μ,s) never equalled the C-bound. — `results/q3_results.json`, Fig. `unlabeled_fig3_information_hierarchy.png`

| weights | true error | sharp \| p + full unlabeled | sharp \| p + pairwise | sharp \| (μ,s) | C-bound | 2μ | 4s |
|---|---|---|---|---|---|---|---|
| uniform | 0.1318 | 0.1844 | 0.2443 | 0.2466 | 0.4232 | 0.5069 | 0.4404 |
| log-odds | 0.1104 | 0.1624 | 0.1871 | 0.3134 | 0.3513 | 0.4321 | 0.3655 |

### Inferences
- For a *fixed finite* ensemble, the C-bound is never sharp given its own inputs: on the grid it was sharp at 0% of points for every finite K. Its optimality statement (Germain Prop. 13) is about an abstract continuous margin variable. Discreteness of the vote count for finite K, and the linear branch μ + E[D] in the region d < R, are both lost. This is a concrete, publishable gap. "min(C, E[L] + E[D])" is a drop-in improvement that is still label-light, since E[D] is unlabeled.
- I did not find the bound L(MV) ≤ E[L] + E[D] stated in the three papers fetched. Germain et al. compare the C-bound only against 2R(G_Q). It may exist elsewhere, so a novelty check is needed (see Gaps).
- For weighted votes, averaged weighted moments lose a lot relative to individual (p, r) (log-odds: 0.313 vs 0.187). That argues for per-classifier partial-identification bounds over averaged second-order PAC-Bayes quantities when K is small.
- These are oracle (population) comparisons only. PAC-Bayes bounds add KL/estimation terms to all quantities, which this study did not model.

### Gaps
- I could not verify against the original papers whether the finite-K closed form (Thm 6) is exactly the Boros–Prékopa (1989) two-moment bound for "at least r of n events". It is very likely a special case.
- I did not search the wider literature (Lacasse et al. 2007; Laviolette et al. 2017; Viallard et al. 2021) for "R(MV) ≤ R(G) + d".
- Empirical (finite-sample, KL-penalised) versions of the new bounds were not implemented.

---

## Q4. Empirical sanity check (scikit-learn 1.9.1)

### Takeaway
On four small binary tasks with five off-the-shelf classifiers, the max–min rule built from (p̂, unlabeled agreement) mostly reduces to the best single classifier (dictator in 67–100% of splits). It matches the best single classifier's test accuracy and beats equal-weight majority and log-odds on 3 of 4 tasks. It never beats the best single classifier meaningfully, because ensembles barely help on these tasks.
- Plug-in sharp intervals are far too narrow. The dominant uncertainty is sampling error in p̂ from the small labeled set, not partial identification.
- A Hoeffding radius on p alone, with the unlabeled law plugged in, gives 100% coverage with widths 0.13–0.26.
- Pure identification widths, computed with population information, are small (0.008–0.115 for majority with the full unlabeled law).

### Cited Findings
**Setup.**
- Tasks: breast_cancer (n=569), digits even-vs-odd (1797), make_classification with redundant/correlated features (6000; 4 informative, 8 redundant, 3 clusters per class, flip_y 0.05, class_sep 0.7), and noisy moons with 4 correlated extra features (4000).
- Splits: 30 random splits of train 40% / labeled validation 10% / unlabeled 25% / test 25%.
- Classifiers: LogisticRegression, GaussianNB, DecisionTree(depth 3), kNN(7), RandomForest(50 trees, depth 5).
- p̂ comes from the labeled split (Laplace-smoothed for log-odds weights). The pairwise agreement and even-part law m̂ come from the unlabeled split. The rule set is all 81 self-dual threshold rules.
- Uncertainty sets:
  - "plug": the smallest common slack making the set feasible.
  - "conf": Hoeffding on p, on agreements, and an L1 radius on m̂, at δ = 0.05.
  - "p-only conf": Hoeffding on p, minimal slack on m̂.
— [q4_empirical.py](S/q4_empirical.py), [q4b_identification_only.py](S/q4b_identification_only.py), [q4c_p_only_conf.py](S/q4c_p_only_conf.py)

**Test accuracy (mean ± sd over 30 splits) [NUM].** — `results/q4_summary.json`, Fig. `unlabeled_fig5_empirical.png`

| method | breast_cancer | digits even/odd | synthetic_corr | moons_noisy |
|---|---|---|---|---|
| best single (by val acc) | 0.9662 ± 0.0169 | 0.9748 ± 0.0076 | 0.7901 ± 0.0140 | 0.8726 ± 0.0095 |
| equal-weight majority | 0.9611 ± 0.0175 | 0.9484 ± 0.0126 | 0.7556 ± 0.0172 | 0.8699 ± 0.0081 |
| log-odds weights | 0.9611 ± 0.0175 | 0.9601 ± 0.0101 | 0.7778 ± 0.0140 | 0.8699 ± 0.0081 |
| max–min, pairwise, plug-in | 0.9643 ± 0.0165 | 0.9714 ± 0.0111 | 0.7906 ± 0.0142 | 0.8727 ± 0.0095 |
| max–min, full unlabeled, plug-in | 0.9653 ± 0.0160 | 0.9715 ± 0.0123 | 0.7906 ± 0.0142 | 0.8724 ± 0.0096 |
| max–min, pairwise, conf | 0.9627 ± 0.0188 | 0.9737 ± 0.0091 | 0.7906 ± 0.0142 | 0.8728 ± 0.0095 |
| max–min, full, conf | 0.9662 ± 0.0175 | 0.9737 ± 0.0091 | 0.7906 ± 0.0142 | 0.8724 ± 0.0096 |
| best rule in hindsight (test) | 0.9786 ± 0.0120 | 0.9765 ± 0.0071 | 0.7945 ± 0.0120 | 0.8774 ± 0.0074 |

- How often the max–min (plug) rule was a dictator: 0.77 / 0.90 / 1.00 / 1.00 (pairwise) and 0.67 / 0.83 / 1.00 / 1.00 (full).
- How often log-odds equalled majority: 1.00 / 0.00 / 0.00 / 1.00.
- Labeled set sizes: 57 / 180 / 600 / 400.

**Interval coverage of the realised test accuracy [NUM].** — `results/q4_summary.json`, `results/q4c_summary.json`
- *Plug-in sets.*
  - Majority: coverage 0.60 / 0.77 / 1.00 / 0.83 (pairwise), mean widths 0.040 / 0.073 / 0.184 / 0.038. With the full unlabeled law: 0.23 / 0.27 / 1.00 / 0.60.
  - Max–min rule: coverage 0.03 / 0.07 / 0.03 / 0.00. Its interval collapses to p̂_i of a dictator, so any miss is pure label-sampling error.
- *Full Hoeffding sets.* Coverage 1.00 everywhere but wide: majority widths 0.32–0.46, max–min-rule widths 0.14–0.24.
- *p-only Hoeffding + plug-in unlabeled law.* Coverage 1.00 on all four tasks for both majority and the max–min rule. Majority widths were 0.258 / 0.197 / 0.240 / 0.186 and max–min widths 0.230 / 0.147 / 0.133 / 0.163, with ε_p = 0.216 / 0.121 / 0.066 / 0.081.

**Identification-only widths (test split as population, exact p and m) [NUM].** — `results/q4b_summary.json`
- Majority interval width: 0.0185 / 0.0573 / 0.1910 / 0.0419 (pairwise) vs 0.0077 / 0.0204 / 0.1147 / 0.0280 (full unlabeled law).
- Test minus lower bound (full): 0.0040 / 0.0115 / 0.0600 / 0.0111.
- With exact information, the max–min rule achieved 0.9769 / 0.9753 / 0.7923 / 0.8764 (pairwise). For comparison: best single 0.9755 / 0.9748 / 0.7923 / 0.8764, majority 0.9611 / 0.9484 / 0.7556 / 0.8699, log-odds 0.9622 / 0.9625 / 0.7770 / 0.8699, and best rule in hindsight 0.9786 / 0.9765 / 0.7945 / 0.8774.

### Inferences
- The check is consistent with the theory. Real classifiers trained on the same features have strongly positively dependent errors, and in that regime the robust rule is "trust the best single model". Majority and log-odds lose 1–3.5 points on digits and synthetic_corr because of correlated weak voters (GaussianNB at 0.647 on digits).
- For a paper, the binding constraint is the labeled sample size. A practical method should propagate the uncertainty in p̂ (Hoeffding/Clopper–Pearson or PAC-Bayes-kl) and treat the unlabeled law as essentially known. This "p-only conf" design achieved full coverage here.
- These tasks do not exercise the regime where robust ensembling beats the best single model (negatively dependent or complementary errors). A convincing empirical section needs diverse-feature or multi-view ensembles, or LLM-judge style ensembles.

### Gaps
- 30 splits per task. Differences between max–min and best single are within one standard error. No statistical tests were run.
- No experiment with deliberately diverse classifiers (different feature subsets or views), where the robust choice would be a non-trivial committee.
- No covariate-shift experiment (labeled source, unlabeled target), which is the main motivating use case.

---

## Q5. Summary of statements and open problems

### Takeaway
The candidate insights are confirmed and sharpened. Unlabeled data identify the even part of π (exactly, when the prior is unknown) and accuracy of self-dual votes lives in the odd part. This yields a clean pair-splitting LP for sharp accuracy intervals, a K=3 closed form for the max–min rule (ΣT_ij < e_min), and a demonstration that standard PAC-Bayes second-order bounds are not sharp for finite ensembles, with an explicit sharp replacement min{C, E[L] + E[D]} as K → ∞.

### Cited Findings
- **Proved here:**
  - Theorem 1: identified set = E(m,p) when the prior is unknown.
  - Theorem 2: identified functionals are those without odd Walsh weight of degree ≥ 3; only dictators have identified accuracy.
  - Theorem 3: TV characterisation when the prior is known.
  - Case (b) identification when η ≠ ½, and the self-dual accuracy identity.
  - The K=3 observable closed form.
  - The "≤ K fractional pairs" structure.
  - The dual certificate.
  - Theorem 4: the K=3 max–min rule and no randomisation gain for K=3.
  - Theorem 5: sharp min{C, 3μ−2s} as K → ∞, i.e. L(MV) ≤ E[L] + E[D] always.
  - Validity of Theorem 6 and of the weighted value-set closed form.
  — the Q1–Q3 sections above; code in `S/`.
- **Numerically verified only:** sharpness of the finite-K and weighted closed forms (max error ≤ 1.3e−7); K=5 max–min ∈ {dictator, triple majority} on 543/543 positive-dependence instances; no gain of non-self-dual predictors over randomised self-dual ones with the prior unknown (0/160 instances with a strict gain); all empirical numbers.
- **Speculative:** the "level-1 Walsh weight / small committee" explanation of sparse robust rules (level-1 weights: dictator 1.0, (3,1,1,1,1) 0.828, 3-majority and (2,1,1,1,0) 0.75, 5-majority 0.703, (2,2,1,1,1) 0.6875; `results/q1b_walsh.json`); that L(MV) ≤ E[L] + E[D] is new; that the finite-K formula coincides with Boros–Prékopa.

### Inferences
Candidate NeurIPS angle: "Label-efficient certification and robust aggregation of classifier ensembles". The ingredients are:
- (i) the even/odd identification theorem;
- (ii) sharp pair-splitting intervals with a finite-sample version (uncertainty in p dominates);
- (iii) max–min / minimax-regret aggregation, with closed forms for K=3 and sparse committees in general;
- (iv) the result that standard second-order PAC-Bayes bounds are not sharp, plus the improved label-light bound min{C, E[L] + E[D]}, and the finite-K discrete version.

The strongest empirical story probably needs ensembles with complementary errors, or distribution shift where unlabeled target data are cheap.

### Gaps
Open problems:
1. A closed form or a combinatorial algorithm for the pair-splitting LP for weighted votes and general K, and the complexity of max–min over all 2-monotone self-dual rules as K grows (the number of rules grows super-exponentially).
2. Characterise when the max–min rule is a dictator, a triple majority or a larger committee (K=5 phase diagram; general K).
3. A finite-sample theory combining labeled uncertainty in p with unlabeled uncertainty in m, including optimal allocation of the labeling budget and PAC-Bayes-kl versions of min{C, E[L] + E[D]} and Theorem 6.
4. Minimax regret vs max–min: closed forms and statistical guarantees.
5. Extra label-free or cheap information (class-conditional prediction rates, a known prior), non-self-dual rules with bias, multiclass (where f_i = f_j no longer equals Z_i = Z_j), and abstention.
6. Literature checks: Boros–Prékopa closed forms; whether R(MV) ≤ R(G) + d appears in the C-bound literature; relation to the Balsubramani–Freund transductive game, which already gives minimax-optimal clipped predictors from unlabeled data but does not treat sharp evaluation of a *given* vote. — [Balsubramani & Freund 2015](https://arxiv.org/pdf/1503.01811)
