# Completion of the thread "partial identification of classifier dependence from labeled marginal accuracies + unlabeled agreement data: sharp bounds and robust weights for weighted votes"

Builds on `unlabeled_agreement_sharp_bounds.md` (earlier session; theorems T1–T6 referenced below are from there) and `ml_bounds_unlabeled_llm.md` (novelty map). Nothing earlier was deleted.

Code folder `S = code/unlabeled_thread/`.
- **New scripts this session:** `q5_identification_extras.py`, `q6_subset_conjecture_and_K5_certs.py`, `q7_finite_sample_and_symK5.py`, `q8_real_data.py`, `q8b_label_regimes.py`, `q9_rate_and_dualnorm.py`, `q10_subset_formula_K7.py`, `q11_label_budget_crossover.py`, `q12_symK5_simplified.py`.
- **Outputs:** `S/results/q5_results.json` … `q12_results.json`, `q8_summary_*.json`, `q8b_summary_*.json`, `q8_records_*.json`, plus `*_log.txt`.
- **Real data:** six PMLB copies of UCI/benchmark tasks (adult, magic, spambase, phoneme, twonorm, churn), downloaded to `code/unlabeled_thread/data_real/`（或环境变量 `PMLB_DATA_DIR`；数据不随仓库分发） from `media.githubusercontent.com/media/EpistasisLab/pmlb/...`.
- **Blocked hosts:** OpenML, UCI archive and HuggingFace are blocked by the session's egress policy (proxy 403). LLM-judge data could therefore not be fetched.
- **Software:** scipy 1.18.1 (HiGHS), numpy 2.5.3, scikit-learn 1.9.1.

Notation is as before. Z_i = 1{classifier i correct}; π = law of Z on {0,1}^K; p_i = P(Z_i=1); e_i = 1−p_i; φ = law of the prediction vector f (from unlabeled data); m(z) = φ(z)+φ(z̄) = π(z)+π(z̄); E(m,p) = {π ≥ 0 : π(z)+π(z̄) = m(z), Zᵀπ = p}. P(p,r) is the set using pairwise joint successes r_ij instead of m. T_ij = P(i and j both wrong), D_ij = P(f_i ≠ f_j). A rule g is self-dual if g(v̄) = 1−g(v), which covers every positive-weight vote without ties. For self-dual g, acc(g) = g·π.

Labels:
- **[PROVED]** means a proof is given here or in the earlier notes.
- **[NUM: script, n, max err]** means numerically verified.
- **[CONJ]** means a conjecture, with the evidence stated.
- **[KNOWN: ref]** means it is already in the literature.

---

## Q1. Exactly what does (labeled marginals + unlabeled prediction law) identify, with and without a class prior, and when do the bounds collapse to a point?

### Takeaway
With the class prior unknown, the identified set for π is exactly the polytope E(m,p). In Walsh/Fourier terms, unlabeled data pin all even-degree correctness moments, the labels pin the degree-1 moments, and the odd moments of degree ≥ 3 are free subject only to positivity. Accuracy of any self-dual vote depends only on odd moments. Hence, with full-support π, no non-dictator vote's accuracy is ever point-identified.

A known prior η adds only two TV constraints. Together with the unlabeled law, η makes the class-conditional accuracies (sensitivity and specificity) identified for free, but the vote-accuracy interval shrinks by only about 5–8%.

Collapse to a point happens only at the boundary. For K=3 the width is exactly 2·(min_{even #errors} π + min_{odd #errors} π), so the interval is a point iff π has a zero on some even-error pattern AND a zero on some odd-error pattern. Even-error patterns are "0 errors" and "exactly two err (3 patterns)"; odd-error patterns are "exactly one errs (3 patterns)" and "all three err". In practice this means some "two-wrong" cell is empty and either the all-wrong cell or some one-wrong cell is empty, as in disjoint-error ensembles. Errors independent of Y with η ≠ ½ also give point identification.

### Cited Findings
**Identification of the set.**
- **Theorem 1 (prior unknown) [PROVED].** The identified set is E(m,p). Pairwise agreement P(f_i=f_j) = P(Z_i=Z_j), so r_ij is identified. — earlier notes, `q1_identification.py`
  - Direct (α,β)-LP vs E(m,p)-LP: max gap ≤ 2e−15 for K = 3, 4, 5. [NUM]
  - The pairwise identity itself is [KNOWN: Platanios, Blum & Mitchell, UAI 2014](https://www.auai.org/uai2014/proceedings/individuals/313.pdf).
  - The same polytope is LP (7) / Prop. 6 of [Mazzetto, Sam, Park, Upfal & Bach, AISTATS 2021](https://par.nsf.gov//servlets/purl/10273522). Its constraint is p_a + p_{r(a)} = h_a + h_{r(a)}, with r = bit flip and h = output law. A fetch of the full text found that the paper asserts the worst case equals the LP but gives no derivation in the main text and does not discuss realisability or point identification. Theorem 1's proof (per-antipodal-pair splitting of α) supplies that sharpness argument.
- **Theorem 2 (what is point-identified) [PROVED; NUM in `q1b_walsh_check.py`, 200+200 instances, err 2.2e−16].**
  - For full-support π, the direction space of E(m,p) is the span of odd Walsh characters of degree ≥ 3.
  - A linear functional c·π is identified iff ĉ has no weight on odd degrees ≥ 3.
  - Among self-dual rules, only the K dictators are point-identified.
  - dim E(m,p) = 2^{K−1} − K; under pairwise-only information the dimension is 2^K − 1 − K − C(K,2). For K=3 the two coincide.
- **Theorem 3 (prior known) [PROVED + NUM, ≤ 1.7e−15].** Add TV(π,φ) ≤ 1−η and TV(π, φ∘flip) ≤ η.

**New this session.**
- **Class-conditional accuracies are identified once η is known [PROVED].**
  - P(f_i=1) = η·a1_i + (1−η)(1−a0_i) and p_i = η·a1_i + (1−η)a0_i, where a1_i and a0_i are the accuracies on Y=1 and Y=0. Two linear equations give (a1_i, a0_i).
  - Consequence: labeled class-conditional accuracies add nothing beyond (p, η, φ). The mean majority-interval widths are exactly equal ("full + η" vs "class-conditional + η + full"): 0.1413 vs 0.1413 for K=3 (n=270) and 0.1328 vs 0.1328 for K=5 (n=125).
  - The interval was never a point (0/270 and 0/125). — `q5_identification_extras.py` (c), `results/q5_results.json`
- **K=3 width formula [PROVED + NUM].**
  - The 1-D direction is δ = χ_{123}(s) = (−1)^{#errors}, and the feasible step range is t ∈ [−min_{δ=+1}π, min_{δ=−1}π].
  - Majority's Walsh/table inner product with δ is 1 − 3 = −2. Hence width(acc_maj) = 2·(min_{#err even} π + min_{#err odd} π).
  - Matched the LP over 3000 instances (4 families including the "anti"/disjoint-error family), max error 5.9e−16. Width was zero in 25.0% of instances, which equals the share of anti-family draws (every 4th instance). The per-family breakdown was not logged separately, but this is consistent with the earlier `q4c` finding of 200/200 zero widths in the anti family. — `q5` (a)
- **Errors independent of Y, η ≠ ½ [PROVED + NUM 1.5e−14].** π = ½[m + (φ−φ̃)/(2η−1)] is fully identified. This is consistent with the classical result that conditionally independent classifiers are identifiable from agreements — [KNOWN: Platanios et al. 2014](https://www.auai.org/uai2014/proceedings/individuals/313.pdf); [Jaffe, Nadler & Kluger, AISTATS 2015](https://ar5iv.arxiv.org/html/1407.7644).

**Value of each information layer (random class-dependent joint laws, p_i > ½) [NUM: `q5` (c)].** "max-min" is the max over all self-dual weighted rules of the worst-case accuracy.

| K (n) | majority width: p only / + pairwise / + full unlabeled / + η | majority lower bound: p only / pw / full / class-cond. | true majority acc | max-min value: p only / pw / full | max-min > best single (pw) |
|---|---|---|---|---|---|
| 3 (270) | 0.418 / 0.149 / 0.149 / 0.141 | 0.570 / 0.717 / 0.717 / 0.719 | 0.785 | 0.771 / 0.777 / 0.777 | 23.3% |
| 5 (125) | 0.475 / 0.263 / 0.144 / 0.133 | 0.525 / 0.693 / 0.755 / 0.760 | 0.829 | 0.798 / 0.817 / 0.818 | 56.8% |

- With p only, the max-min value equals the best single p_max in every instance. This is consistent with [KNOWN: Mazzetto et al. 2021 Prop. 3, "given only error rates, no function beats the best single labeler"](https://par.nsf.gov//servlets/purl/10273522) and with [de Oliveira, Ishii & Lin, arXiv 2106.00088](https://arxiv.org/abs/2106.00088): binary robust aggregation relies on a single source.

### Inferences
- The clean statement for a paper is the following. *(φ, p) identifies exactly the even part plus degree-1 moments; every vote's accuracy is a known affine term plus an unidentified odd-degree-≥3 term, bounded only by positivity.* "Known η" and "class-conditional accuracies" are essentially worthless additions in case (a): they cut width by 5–8%.
- Point identification is a measure-zero, boundary event, except under structural assumptions (independence from Y with η ≠ ½). The anti family's zero widths come from empty cells, not from information.

### Gaps
- For general K, a full characterisation of point identification at boundary π (the face of E containing π lying in a level set of g) was not turned into a combinatorial criterion beyond K=3.
- Multiclass labels were not analysed. There, P(f_i=f_j) ≠ P(Z_i=Z_j), because two classifiers can be wrong in different ways.

---

## Q2. Closed forms (K=3, K=5), max-min weights, and size of improvement over marginal-only and over C-bound / second-order bounds

### Takeaway
**K=3.** Everything is closed form. The sharp majority error is min{ΣT_ij, min_i(e_i+d_i)}, and the max-min rule is "majority iff Σ_{i<j}T_ij < e_min". Both are algebraically identical to Mazzetto et al. 2021 Prop. 5, so they are **known**.

**K=5, symmetric (p_i = p, T_ij = T).** Fully closed. Max-min error = **min{e, 3T, 1 − 5e + 10T}** (dictator / any triple / full majority).

**K=5, general pairwise information.** The max-min value equals **max{p_max, max over the 10 triples of the K=3 closed form, L(5-majority)}** to 1e−15. Only L(5-majority) still needs a 32-variable LP: harvesting dual vertices produced 1849 certificates and still did not give a closed form (exact in only 33.6% of fresh instances).

**Weighted vs subset votes.** A central new observation: under pairwise information, no weighted vote ever beat the best odd-subset majority.
- Checked with adversarial search for K=4 and K=5 and random plus local search for K=7.
- So "robust weights" collapse to robust subset selection, which is Mazzetto's algorithm.
- With the full unlabeled law, weighted rules can strictly win (adversarial gain up to 0.067 for K=4 and 0.015 for K=5), but never did on natural families.

**Comparison with existing bounds.**
- Synthetic (earlier notes): the sharp bound beats min{FO, TND, C-bound} by a mean of 0.17–0.22 for K=3 and 5.
- Real data (this session, oracle population quantities, K=5): mean majority-error bound 0.202 (sharp | p + pairwise) vs C-bound 0.446, FO 0.377, TND 0.481, true 0.150.

### Cited Findings
**K=3 closed forms.**
- **Sharp interval [PROVED + NUM].** A ∈ [max{max_i(p_i−d_i), Σp−2+D}, min{min_i(p_i+d_i), Σp−D}], with d_i = (D_ij + D_ik − D_jk)/2. Max-min picks majority iff ΣT < e_min. — earlier notes; `q3b_closed_forms.py` (2000/2000)
- **Equivalence to Mazzetto Prop. 5 [NUM: `q5` (b), 3000 instances incl. anti family, max err 5.0e−16].**
  - Prop. 5 states worst-case error = Σe − (D12+D13+D23)/2 − 2q, with q = ½max{e_i+e_j−D_ik−D_jk (three choices), 0}.
  - It equals our min{ΣT, min_i(e_i+d_i)} and the LP.
  - The symmetric threshold (majority robust iff d ≥ 4ε/3) is stated in the same paper; it is equivalent to our ρ*(p) threshold, since 3T < e ⟺ d > 4e/3 because T = e − d/2.
  - [KNOWN: Mazzetto et al. AISTATS 2021](https://par.nsf.gov//servlets/purl/10273522)

**K=5, symmetric closed form [PROVED for the rule values + NUM that no other type wins].**
- **5-majority (from T6 with Eν = 5e, E[ν(ν−1)] = 20T).** Worst error = min{1, (35e−20T)/15, 10T/3, 1−5e+10T}.
- **Triple:** 3T. **Dictator:** e.
- (35e−20T)/15 ≥ e and 10T/3 ≥ 3T always hold, so the max-min error is min{e, 3T, 1−5e+10T}.
  - Checked against the full LP grid for all 7 rule types: 11,680 feasible (p,ρ) cells, max error 7.5e−15. The argmin was the dictator in 9029 cells, the triple in 2091 and the 5-majority in 560. — `q12_symK5_simplified.py`, `results/q12_results.json`
  - The LP itself matched T6 (max err 1.1e−15) and the K=3 form (2.2e−16) on 858 cells. — `q7` (B)
- **Phase boundaries [PROVED from the formula].**
  - The triple beats the dictator iff T < e/3. Under independence (T = e²) this means p > 2/3.
  - The 5-majority wins iff T < min{(5e−1)/7, (6e−1)/10}. This needs e > 1/5, and it is impossible under independence because 1−5e+7e² > 0 always.
  - **So under exact independence the max-min rule never uses all 5 voters.**

**K=5, general pairwise information.**
- Max-min value = max{p_max, max_triples A3, L5(maj)}. 300 instances (incl. anti and sparse-Dirichlet laws), max error 1.1e−15. — `q10_subset_formula_K7.py` (1)
- The worst case of a triple majority inside the 5-voter set equals the 3-voter closed form on that triple: 3760 triple-instances, max |diff| 1.7e−15. So sub-committees "project" exactly. — `q5` (d)
- **No compact closed form for L5(maj) [NUM, negative].**
  - Harvesting optimal dual vertices of the L5 LP from 1500 instances gave 40 S₅-orbits (1849 certificates).
  - "max over harvested certificates" was exact in only 336/1000 fresh instances, max gap 0.014.
  - The vertex set is large, so the LP (32 variables) is the practical tool. — `q6_subset_conjecture_and_K5_certs.py`, `results/q6_results.json`

**Weighted vs subset votes.** Search maximised gain = maxmin(all self-dual LTFs) − maxmin(odd-subset majorities). — `q6`, `q10`, `q5` (e)

| setting | info | search | instances (LPs) | max gain |
|---|---|---|---|---|
| K=4 (12 rules, 8 subset) | pairwise | adversarial hill-climb over π | 2460 | **0** |
| K=4 | full unlabeled | same | 2460 | **0.0667** |
| K=5 (81 rules, 16 subset) | pairwise | adversarial | 525 | **0** |
| K=5 | full unlabeled | adversarial | 525 | 0.0152 |
| K=5, 5 natural families (indep., mixture, bloc, anti, exchangeable) | pairwise and full | — | 376 | 0 (≤ 1e−15) |
| K=7 (64 subset rules vs 150 random + 3×60-step local-search weights) | pairwise | random + local search | 40 | 0 (2e−14) |
| K=7 | full unlabeled | same | 25 | 0 |

**Comparison with existing bounds.**
- **Synthetic, from the earlier notes [NUM, `q3_ml_bounds.py`].** Over the (μ,s) grid, the best existing bound min{2μ, 4s, C} exceeds the sharp bound by a mean of 0.222 (K=3), 0.168 (K=5) and 0.095 (K→∞). The C-bound is never sharp for finite K. The sharp continuum bound is min{C, E[L]+E[D]}.
- **Real data, oracle population quantities on the test split [NUM, `q8_real_data.py`].** 6 datasets × 2 ensembles × 10 splits; uniform 5-majority.
  - "same" = 5 sklearn models on all features. "views" = each model on its own disjoint random feature block.

| task/ensemble | true err | FO 2E[L] | TND 4E[LL′] | C-bound | min(C, E[L]+E[D]) | T6 sharp \| (μ,s) | sharp \| p only | sharp \| p + pairwise | sharp \| p + full |
|---|---|---|---|---|---|---|---|---|---|
| adult/same | 0.157 | 0.350 | 0.520 | 0.485 | 0.265 | 0.250 | 0.292 | 0.180 | 0.171 |
| adult/views | 0.222 | 0.454 | 0.659 | 0.600 | 0.352 | 0.331 | 0.379 | 0.270 | 0.257 |
| magic/same | 0.168 | 0.409 | 0.557 | 0.527 | 0.335 | 0.313 | 0.341 | 0.217 | 0.200 |
| magic/views | 0.270 | 0.571 | 0.835 | 0.734 | 0.439 | 0.413 | 0.476 | 0.349 | 0.319 |
| twonorm/same | 0.023 | 0.150 | 0.114 | 0.112 | 0.112 | 0.056 | 0.058 | 0.031 | 0.027 |
| twonorm/views | 0.060 | 0.410 | 0.297 | 0.271 | 0.271 | 0.139 | 0.341 | 0.139 | 0.077 |
| spambase/same | 0.071 | 0.231 | 0.235 | 0.235 | 0.226 | 0.148 | 0.192 | 0.121 | 0.110 |
| spambase/views | 0.112 | 0.420 | 0.374 | 0.366 | 0.358 | 0.215 | 0.342 | 0.213 | 0.185 |
| phoneme/same | 0.199 | 0.421 | 0.594 | 0.554 | 0.335 | 0.314 | 0.351 | 0.246 | 0.226 |
| phoneme/views | 0.271 | 0.590 | 0.771 | 0.709 | 0.499 | 0.442 | 0.491 | 0.374 | 0.349 |
| churn/same | 0.107 | 0.236 | 0.342 | 0.329 | 0.183 | 0.172 | 0.197 | 0.135 | 0.125 |
| churn/views | 0.141 | 0.278 | 0.479 | 0.434 | 0.178 | 0.171 | 0.232 | 0.144 | 0.143 |
| **mean** | **0.150** | 0.377 | 0.481 | 0.446 | 0.296 | 0.247 | 0.308 | **0.202** | **0.182** |

- These are oracle forms. PAC-Bayes versions add KL and estimation terms to every column, and the earlier session did not model these.
- The ordering sharp|full ≤ sharp|pw ≤ T6 ≤ min(C, E[L]+E[D]) ≤ C holds in every row.
- FO, TND and the C-bound follow [Wu, Masegosa, Lorenzen, Igel & Seldin, NeurIPS 2021](https://arxiv.org/pdf/2106.13624.pdf), [Masegosa et al., NeurIPS 2020](https://ar5iv.labs.arxiv.org/html/2007.13532) and [Germain et al., JMLR 2015](https://ar5iv.arxiv.org/html/1503.08329).

### Inferences
- **Improvement over marginal-only.** For K=3, pairwise or unlabeled information cuts the majority interval from about 0.42 to 0.15 wide. For K=5 it cuts from 0.47 to 0.26 (pairwise) and 0.14 (full law). On real data the error bound falls from 0.308 (p only) to 0.202 (pairwise) and 0.182 (full).
- **Improvement over second-order PAC-Bayes oracles.** About a 2.2× smaller bound than the C-bound on average on real data. Even T6, which uses only averaged (μ,s) but exploits finite-K discreteness, beats the C-bound everywhere.
- **The big structural message is negative for "weights".** Given (p, pairwise agreement), robust aggregation = choose the best odd sub-committee and take its plain majority, empirically for K ≤ 7. If proved, this says weighting is useless under pairwise ambiguity. That is a clean theorem, but it makes the "robust weights" headline essentially Mazzetto's subset-selection procedure. Weights matter only when higher-order unlabeled information is used, and even then rarely.

### Gaps
- No proof of the subset-majority conjecture (see Q6, theorem candidate C). An LP-duality or "pairwise-polytope extreme points" argument is the likely route.
- No closed form for asymmetric L5(maj). Mazzetto et al. also give only the K=3 closed form.
- I did not verify against the original that T6 coincides with the Boros–Prékopa (1989) closed form. The publisher page only says the paper gives "closed form lower and upper bounds… when the basis is dual feasible" for at-least-r probabilities via binomial moments ([Boros & Prékopa, Math. Oper. Res. 14(2):317, 1989](https://pubsonline.informs.org/doi/fpi/10.1287/moor.14.2.317)). The r=1 S₁/S₂ sharp bounds (Kwerel, Dawson–Sankoff) are summarised in [Boros, RUTCOR notes](https://rutcor.rutgers.edu/Boros.pdf).

---

## Q3. Finite-sample version: consistency, confidence-set LP, rate — is a clean theorem within reach?

### Takeaway
Yes, and it is short.
- Define L̂(g) as the minimum of g·π over the confidence set C = {π ∈ Δ : ‖Zᵀπ − p̂‖_∞ ≤ ε_p, Σ_{pairs}|π(v)+π(v̄) − m̂(v)| ≤ ε_m}.
- With probability ≥ 1−δ: **L(g) − 2(ε_p‖λ*_g‖₁ + ε_m‖μ*_g‖_∞) ≤ L̂(g) ≤ L(g) ≤ acc(g)**.
  - (λ*, μ*) is any optimal dual of the population LP.
  - Its norms are bounded by a constant C_K(g) that depends only on (K, g), because the constraint matrix depends only on K.
  - ε_p = √(log(4K/δ)/(2n_ℓ)) and ε_m = √(2(2^{K−1} log 2 + log(4/δ))/n_u).
- So the certificate is valid and converges at rate n^{−1/2}. The max-min rule chosen on C has population worst case within 2·max_g C_K(g)(ε_p+ε_m) of the optimum.
- Numerically, validity and the gap inequality held in 100% of cases, and the empirical log–log slope was −0.44.
- **Practical insight ("price of marginal labels").** A dictator's certificate loses ε_p (‖λ‖₁ = 1). A K-committee's certificate loses up to about ‖λ‖₁·ε_p, with ‖λ‖₁ up to K for majority. So with marginal-only labels a committee is chosen only if its identified advantage over the best single exceeds about (‖λ‖₁ − 1)ε_p.
- On real data this kept the confidence-set max-min rule a dictator in 120/120 split-ensemble cases at n_ℓ = 200 (8 to 10 splits for each of the 12 settings).

### Cited Findings
**Proof sketch of the theorem [PROVED].**
1. Write the population LP as V(b) = min{g·π : Aπ = b, π ≥ 0}, with b = (m, p). By LP duality V(b) = max{b·y : Aᵀy ≤ g}. This is a maximum of linear functions of b, so it is convex, and any dual optimum y* is a subgradient: V(b′) ≥ V(b) + y*·(b′−b).
2. On the event ‖p̂−p‖_∞ ≤ ε_p and ‖m̂−m‖₁ ≤ ε_m (Hoeffding + union bound, and the Weissman et al. L₁ multinomial bound on the 2^{K−1}-cell law m), true π ∈ C, so L̂ ≤ L.
3. Any π ∈ C has b′ = Aπ with ‖b′_p − b_p‖_∞ ≤ 2ε_p and ‖b′_m − b_m‖₁ ≤ 2ε_m. Hence g·π ≥ V(b′) ≥ V(b) − 2ε_p‖λ*‖₁ − 2ε_m‖μ*‖_∞ (Hölder).
4. A has full row rank 2^{K−1}+K, since dim E = 2^K − rank = 2^{K−1} − K. So the dual polyhedron is pointed, y* can be taken at one of finitely many vertices determined by (K, g), and C_K(g) := max over vertex dual norms is finite.
5. Selection: L(ĝ) ≥ L̂(ĝ) ≥ L̂(g*) ≥ L(g*) − 2max_g C_K(g)(ε_p+ε_m). ∎

**Numerical checks [NUM].**
- *Coverage and gap, K=5, 38 random joint laws × 3 sample sizes × 81 rules.* `q7_finite_sample_and_symK5.py` (A), `results/q7_results.json`. In every case the event held, L̂ ≤ L held, and the gap was ≤ the bound (100%).

| (n_ℓ, n_u) | ε_p | ε_m | mean gap L−L̂ (all rules) | majority gap | mean theoretical bound | mean regret of chosen max-min rule (max) |
|---|---|---|---|---|---|---|
| (200, 5000) | 0.122 | 0.079 | 0.243 | 0.282 | 1.025 | 0.029 (0.097) |
| (800, 20000) | 0.061 | 0.039 | 0.145 | 0.177 | 0.513 | 0.023 (0.099) |
| (3200, 80000) | 0.031 | 0.020 | 0.085 | 0.112 | 0.256 | 0.022 (0.094) |

- *Rate.* Majority gap with n_ℓ = 3200 → 204,800 and n_u = 25n_ℓ (12 laws): 0.110, 0.070, 0.035, 0.018. The log–log slope was −0.44, close to −½. Validity and the bound held in all cases. — `q9_rate_and_dualnorm.py`, `results/q9_results.json`
- *Dual-norm constants (observed maxima of ‖λ*‖₁ / ‖μ*‖_∞).* These are not proven constants.
  - All rules: K=3: 3 / 2; K=4: 5 / 4; K=5: 7 / 6; K=6: 5 / 4; K=7: 7 / 6.
  - Majority: ‖λ*‖₁ ≤ K observed (3, 5, 7).
  - By type for K=5 (q7): dictator 1/0, triple 3/2, every other type ≤ 5/4. — `q9` (900, 750, 900, 300, 360 LPs)

**Relation to econometric inference on LP-defined bounds.**
- Value functions of LPs are only directionally differentiable, so the naive bootstrap can fail. [Hsieh, Shi & Shum (working paper)](https://users.ssc.wisc.edu/~xshi/research/lprog.pdf) note this, citing Fang & Santos (2019), and give KKT-based uniformly valid confidence sets for LP solutions.
- [Fang, Santos, Shaikh & Torgovitsky, Econometrica 91(1) 2023](https://ideas.repec.org/p/arx/papers/2009.08568.html) test whether a large linear system with known coefficients has a non-negative solution. That test applies directly to "is E(m̂,p̂) consistent", i.e. a specification test that the reported accuracies are compatible with the unlabeled law.
- Our theorem is the simple, conservative, finite-sample (non-asymptotic) counterpart. It does not achieve the sharper asymptotic size control of that literature.

**Plug-in consistency.**
- V is Lipschitz on its feasible domain, so L(p̂, m̂) → L(p, m) whenever the plug-in set is non-empty.
- The plug-in set can be empty with positive probability when the true (p, m) lies on the boundary of the feasible region, as with the anti family. The earlier session handled this with a "minimal slack".
- Plug-in intervals had poor coverage: 0–83% for the max-min rule's interval in the earlier `q4`. Only the confidence-set version is a valid theorem.

### Inferences
- The clean theorem is within reach and essentially complete: valid finite-sample certificates and O(√(log(K/δ)/n_ℓ) + √(2^{K−1}/n_u)) suboptimality of the max-min rule.
- The 2^{K−1} term means the full-law version needs n_u ≫ 2^K; the pairwise version replaces 2^{K−1} by C(K,2)·log terms. This is exactly the classical efficiency/robustness trade-off, and it argues for the pairwise set when K ≳ 8.
- The ‖λ*‖₁ factor is the key interpretable quantity. It explains why the earlier `q4` and the new `q8` experiments almost always return a dictator at realistic label budgets. A paper should state it as a proposition: a committee g is certified-preferred over the best single only if L(g) − p_max > (‖λ*_g‖₁ − 1)ε_p, up to ε_m terms.
- An "optimal allocation of the labeling budget" result follows directly: minimise Σ_i |λ*_i|ε_{p,i}, which puts labels on classifiers with large dual weight.

### Gaps
- No proof of a closed-form bound on C_K(g), for example C_K(maj) ≤ K, which matches every observed value.
- No lower bound or minimax rate showing that n^{−1/2} is optimal.
- Not derived: PAC-Bayes-kl or empirical-Bernstein radii (tighter than Hoeffding when p ≈ 1); confidence sets for the pairwise version (Hoeffding on C(K,2) agreements was used in `q4`).

---

## Q4. Empirical evaluation: what did q4 do, and what does a real experiment show?

### Takeaway
The earlier `q4_empirical.py` used two small real sklearn-bundled datasets (breast_cancer, digits even/odd) and two synthetic ones (make_classification, noisy moons) with 5 sklearn classifiers. It was a sanity check.

This session ran a real-data experiment on six PMLB/UCI tasks (adult 48,842 rows; magic 19,020; spambase 4,601; phoneme 5,404; twonorm 7,400; churn 5,000). It used two K=5 ensembles ("same features" and "disjoint feature views") and 10 splits each. Findings:
- (i) Confidence-set certificates were valid for all 81 rules in 120/120 split-ensembles.
- (ii) Sharp bounds are far tighter than C-bound, TND and FO (table in Q2).
- (iii) With 200 marginal labels, the max-min rule was **always the best single classifier**. It therefore loses to majority or log-odds exactly where ensembling helps (twonorm/views: 0.818 vs 0.940; spambase/views: 0.858 vs 0.889) and wins where weak correlated voters hurt (phoneme/same: 0.832 vs 0.801).
- (iv) Even with exact population information, the max-min rule was a dictator in 11/12 settings. A committee became robust-optimal only on twonorm/views, and only once ε_p ≤ 0.03, about 3,300 labels per classifier.
- (v) If labels are **joint** (one shared labeled set), simply estimating every rule's accuracy with a union bound does as well or better than partial identification.

So the method's practical niche is the **marginal-only-labels** regime: separate gold sets per annotator or model, vendor-reported accuracies, or privacy-separated evaluation.

### Cited Findings
**Setup [`q8_real_data.py`].**
- n_ℓ = 200 marginal labels; n_u = min(5000, 0.4·rest) unlabeled; train ≤ 4000; the rest is test.
- Models: LR, GaussianNB, Tree(d=3), kNN(7), RF(50, d=5).
- Rule set: all 81 self-dual LTFs.
- Confidence set: Hoeffding ε_p = 0.122; Weissman ε_m ≈ 0.077 at n_u = 5000.

Test accuracy (mean over 10 splits):

| task/ensemble | best single | majority | log-odds | max-min (plug-in) | max-min (conf) | best rule in hindsight | certified LB of max-min (conf) | certified LB of majority |
|---|---|---|---|---|---|---|---|---|
| adult/same | 0.8426 | 0.8433 | 0.8433 | 0.8405 | 0.8475 | 0.8483 | 0.724 | 0.666 |
| adult/views | 0.8039 | 0.7778 | 0.7778 | 0.8039 | 0.8039 | 0.8128 | 0.691 | 0.592 |
| magic/same | 0.8387 | 0.8315 | 0.8316 | 0.8387 | 0.8387 | 0.8414 | 0.730 | 0.649 |
| magic/views | 0.7578 | 0.7299 | 0.7366 | 0.7578 | 0.7578 | 0.7665 | 0.647 | 0.520 |
| twonorm/same | 0.9766 | 0.9767 | 0.9776 | 0.9762 | 0.9766 | 0.9797 | 0.859 | 0.810 |
| twonorm/views | 0.8174 | **0.9401** | 0.9401 | 0.9076 | 0.8181 | 0.9401 | 0.714 | 0.567 |
| spambase/same | 0.9185 | 0.9291 | 0.9291 | 0.9209 | 0.9185 | 0.9319 | 0.804 | 0.701 |
| spambase/views | 0.8581 | 0.8885 | 0.8880 | 0.8772 | 0.8581 | 0.8997 | 0.745 | 0.579 |
| phoneme/same | **0.8315** | 0.8008 | 0.8076 | 0.8315 | 0.8315 | 0.8339 | 0.722 | 0.591 |
| phoneme/views | 0.7411 | 0.7287 | 0.7319 | 0.7373 | 0.7373 | 0.7623 | 0.627 | 0.460 |
| churn/same | 0.9031 | 0.8929 | 0.8929 | 0.9000 | 0.9001 | 0.9044 | 0.783 | 0.698 |
| churn/views | 0.8713 | 0.8592 | 0.8592 | 0.8713 | 0.8701 | 0.8754 | 0.760 | 0.694 |

- Coverage of realised test accuracy by the confidence interval, averaged over all 81 rules: 1.000 in every setting.
- The confidence-set max-min rule was a dictator in 10/10 splits in every setting.
- **Note on the "max-min (conf)" column.** It is not always the "best single" column: max-min chooses the dictator with the highest certified lower bound, which can differ from the validation argmax. For example, adult/same 0.8475 vs 0.8426.
- The plug-in max-min rule sometimes picks a committee (twonorm/views 0.908; spambase/views 0.877). Its intervals are not valid (earlier `q4`: 0–7% coverage).
- The best rule in hindsight is a non-dictator in most "views" settings, for example the 5-majority 10/10 on twonorm/views.

**Labeled-information regimes [`q8b_label_regimes.py`, 8 splits].**
- MARG = marginal labels (each classifier's accuracy on its own n_ℓ points) + unlabeled law → confidence-set max-min.
- JOINT = one shared labeled set of n_ℓ points, all 81 rule accuracies estimated directly, Hoeffding + union-bound certificate.
- HYBRID = JOINT constraints + unlabeled law in one LP.

Mean over 12 settings at n_ℓ = 200:

| regime | test acc | certified LB | certificate coverage | rule types chosen |
|---|---|---|---|---|
| oracle | 0.8617 | — | — | — |
| MARG | 0.8352 | 0.7340 | 1.00 | dictator 96/96 |
| JOINT | 0.8541 | 0.7313 | 1.00 | mixed; committees in most settings |
| HYBRID | 0.8517 | 0.7413 | 1.00 | mostly dictator, committees on "views" |

- Biggest gaps: twonorm/views MARG 0.808 vs JOINT 0.925; spambase/views 0.824 vs 0.871.
- At n_ℓ = 1000 (4 large-data settings): MARG 0.8137, JOINT 0.8141, HYBRID 0.8146, oracle 0.8167. Certified LBs: MARG 0.764 > HYBRID 0.758 > JOINT 0.752.

**Label budget needed for a committee [`q11_label_budget_crossover.py`].** Test set used as population; unlabeled law exact; ε_p swept from 0 to 0.137.
- 11 of 12 settings: the dictator is max-min optimal even at ε_p = 0.
- twonorm/views: 5-majority at ε_p ≤ 0.01 (acc 0.9368, cert 0.869); triple at ε_p = 0.03 (acc 0.902); dictator at ε_p ≥ 0.06 (acc 0.809). The ε_p = 0.03 point corresponds to n_ℓ ≈ log(400)/(2·0.03²) ≈ 3,300 labels per classifier.

### Inferences
- The method is valid but conservative. On off-the-shelf tabular ensembles its robust choice is "the best single model", and it is right to be cautious only when ensembling actually hurts (phoneme/same, magic, adult/views).
- To demonstrate value, an experiment must sit in a regime where (a) only per-classifier accuracies are available and (b) the identified set is narrow enough. That means many marginal labels per voter (≥ 10³), a full unlabeled law, and complementary voters.
  - **Crowdsourcing with per-worker gold questions** fits naturally: worker accuracy comes from that worker's own gold items, and agreements come from overlapping unlabeled items.
  - **LLM judges with published per-judge benchmark accuracies** also fit, but the transfer of p_i across distributions is then an extra assumption.
- Reporting the certified lower bound, not just accuracy, is where the method shines: bounds of 0.63–0.86, always valid, versus PAC-Bayes oracle bounds that are often near-vacuous (Q2 table).

### Gaps
- No LLM-judge or WRENCH experiment was possible: huggingface.co, openml.org and archive.ics.uci.edu are blocked by the session's egress policy.
- Not run: crowdsourcing data with gold questions; covariate shift (labels on source, unlabeled target); B&F HedgeClipper, Dawid–Skene or Snorkel baselines.
- 8–10 splits per setting, no significance tests.

**Proposed experiment design for the team.**
- **E1 (validity).** Coverage and width of certified bounds vs n_ℓ and n_u on UCI plus CIFAR-10 binary deep ensembles. Baselines: PAC-Bayes FO/TND/C-bound with KL terms.
- **E2 (marginal-label regime).** Crowdsourcing datasets with gold questions (e.g., RTE/TEMP/Bluebirds-style binary crowd data), where each worker's accuracy comes from their own gold items. Compare MARG max-min vs MV, Dawid–Skene, Mazzetto subset selection and HedgeClipper.
- **E3 (LLM judges).** 5–7 judges on binary-ised RewardBench/MT-Bench-human pairs, with p_i from a held-out labeled slice per judge and the unlabeled law from the rest. Add a same-family redundant judge and show that certified bounds and the chosen sub-committee react to the redundancy.

---

## Q5. Novelty: what is already known (weak supervision, PAC-Bayes, partial identification)?

### Takeaway
More of this thread is known than the earlier notes assumed.
- **Known:**
  - Mazzetto et al. (AISTATS 2021): the LP over E(m,p) (their LP 7 / Prop. 6) and the pairwise LP (Prop. 4); the K=3 closed form (Prop. 5, verified identical); the marginal-only "no rule beats the best single labeler" (Prop. 3); and robust subset-majority selection (exhaustive or greedy).
  - The pairwise identity and identifiability under independence (Platanios 2014; Jaffe et al. 2015).
  - LP partial identification under measurement error with multiple proxies of unknown dependence: Molinari 2008; Finkelstein et al., UAI 2021.
  - Sharp Fréchet bounds for weak-supervision metrics: Maia Polo et al., NeurIPS 2024.
- **Apparently new (none found in the papers checked):**
  - the Walsh even/odd identification theorem and "only dictators are identified";
  - the sharpness and realisability proof for E(m,p), the TV characterisation with known η, and the class-conditional corollary;
  - the empirical finding that weighted votes never beat subset majorities under pairwise ambiguity;
  - the symmetric K=5 closed form and phase diagram;
  - the sharp min{C-bound, E[L]+E[D]} continuum bound and the finite-K T6 bound as replacements for second-order PAC-Bayes oracles;
  - the finite-sample dual-norm certificate theorem with the "price of marginal labels".

### Cited Findings
**Weak supervision / unsupervised ensembles.**
- [Mazzetto, Sam, Park, Upfal & Bach, "Semi-Supervised Aggregation of Dependent Weak Supervision Sources With Performance Guarantees", AISTATS 2021](https://par.nsf.gov//servlets/purl/10273522), per full-text fetch:
  - Inputs: ε_i from labeled data; D_ij or the output law Q from unlabeled data. "We make no assumptions on the choice of the labelers… and possible dependencies".
  - Equality-constrained LPs over correctness patterns; the K=3 closed form; only unweighted subset majorities (with greedy subset growth).
  - No finite-sample theory, no class-prior discussion, and no point-identification discussion.
- Balsubramani & Freund (COLT 2015) give a transductive minimax over the full unlabeled prediction matrix with accuracy lower bounds and no independence assumption. The optimal predictor is a clipped, randomised weighted vote. — [arXiv 1503.01811](https://arxiv.org/pdf/1503.01811)
- Corrada-Emmanuel, ["No-Knowledge Alarms for Misaligned LLMs-as-Judges"](https://arxiv.org/html/2509.08593):
  - an LP over integer response counts consistent with observed agreements;
  - feasibility "alarms" only: "cannot tell us which of the graders is the one that is misaligned";
  - multiclass, two graders worked out;
  - no vote-accuracy bounds and no finite-sample treatment.
- [Haber & Wexler, "Practical estimation of ensemble accuracy", NeurIPS 2023 InfoCog workshop](https://neurips.cc/virtual/2023/77161): a combinatorial unlabeled bound on ensemble mistakes; partial identification is not mentioned in the abstract.

**PAC-Bayes majority vote.**
- The C-bound is sharp only for a continuous margin variable with μ₂ ≤ μ₁ (Germain et al., Prop. 13) — [JMLR 2015](https://ar5iv.arxiv.org/html/1503.08329).
- [Viallard, Germain, Habrard & Morvant (arXiv 2104.13626)](https://arxiv.org/pdf/2104.13626) discuss 2r = 2e + d and 4e but, per fetch, state no linear "Gibbs + disagreement" bound and no finite-voter sharpness analysis.
- I found no source for L(MV) ≤ E[L] + E[D]; the novelty of that bound is unconfirmed, not established.

**Partial identification (econometrics / causal ML).**
- Molinari, "Partial identification of probability distributions with misclassified data", J. Econometrics 144(1):81–117 (2008) — [IDEAS](https://ideas.repec.org/a/eee/econom/v144y2008i1p81-117.html). This is the general Manski-style framework for misclassification. I could not fetch its content.
- [Finkelstein, Adams, Saria & Shpitser, "Partial identifiability in discrete data with measurement error", UAI 2021](https://ar5iv.labs.arxiv.org/html/2012.12449): a latent variable with multiple proxies whose dependence is left unrestricted, LP bounds sharp when assumptions are linear. Our problem is a special instance: latent Y, proxies f_i, "accuracy" assumptions.
- Inference: [Fang, Santos, Shaikh & Torgovitsky, Econometrica 2023](https://ideas.repec.org/p/arx/papers/2009.08568.html); [Hsieh, Shi & Shum](https://users.ssc.wisc.edu/~xshi/research/lprog.pdf).
- Weak-supervision evaluation via partial identification: [Maia Polo et al., NeurIPS 2024](https://proceedings.neurips.cc/paper_files/paper/2024/hash/f4c6bec746b0aeca8c2cd15096f1ad1f-Abstract.html).

**Novelty map of this thread** (status of each item as reported in this session's notes):

| claim | status | prior work |
|---|---|---|
| r_ij from agreement | KNOWN | Platanios 2014 |
| E(m,p) as LP; pairwise LP | KNOWN (as LP) | Mazzetto 2021 LP 5/7 |
| sharpness / realisability of E(m,p) (T1), TV form with η (T3) | PROVED here; Mazzetto assert worst-case = LP but give no proof in the main text | — |
| even/odd Walsh identification (T2), "only dictators identified" | PROVED; not found elsewhere | sign-ambiguity folklore (Jaffe 2015, Parisi 2014) |
| K=3 closed form + max-min rule | KNOWN | Mazzetto Prop. 5 (verified identical) |
| marginal-only ⇒ dictator | KNOWN | Mazzetto Prop. 3; de Oliveira et al. |
| symmetric K=5 closed form / phase diagram | PROVED (values) + NUM; new | — |
| weighted votes never beat subset majorities under pairwise info | CONJ (strong NUM); new; would make Mazzetto's subset search max-min-optimal over all weighted votes | — |
| full unlabeled law ⇒ weights can help | NUM counterexamples; new | B&F randomised version |
| min{C, E[L]+E[D]} continuum sharp; T6 finite K | PROVED; ML-new as far as checked; T6 likely a special case of Prékopa discrete moment bounds | Germain 2015; Boros–Prékopa 1989 |
| finite-sample certificate + dual-norm "price of marginal labels" | PROVED; simple; new in this ML context | generic LP sensitivity; econometric inference literature |

### Inferences
- A NeurIPS reviewer who knows Mazzetto 2021 will see the E(m,p) LP, the K=3 results and subset-MV selection as prior art. The paper must lead with what Mazzetto lacks:
  - (1) the identification theory (Walsh, only dictators identified, sharpness proof, η and class-conditional corollaries);
  - (2) the weighted-vs-subset theorem;
  - (3) finite-sample certificates with the ‖λ‖₁ price;
  - (4) the sharp replacement for PAC-Bayes second-order oracles.
- The econometric framing (Manski/Molinari, Finkelstein et al.) should be cited as the umbrella. The contribution is a crisp, ML-specific structure theorem plus statistics, not the idea of LP bounds.

### Gaps
- I did not read in full: Molinari 2008; the Mazzetto appendix (to see whether Prop. 6's proof contains the realisability argument); Corrada-Emmanuel's other NTQR papers (2312.05392, 2510.00821); Jaffe et al. 2016 (dependent classifiers).
- No OpenReview search of the 2025–26 cycles was done.

---

## Q6. Candidate NeurIPS-style theorem statements (with novelty risk and 6-week feasibility for 3–4 students)

### Takeaway
Four candidates. The best-balanced package is **A + C + D**:
- identification structure;
- the "weights don't help under pairwise ambiguity" theorem, a conjecture with strong evidence;
- finite-sample certificates with the price of marginal labels.

B, the sharp replacement for PAC-Bayes oracles, is a solid stand-alone or secondary contribution.

### Cited Findings
**Theorem A (Identification structure) — status PROVED.**
- *Statement.* Let K binary classifiers have marginal accuracies p and let the unlabeled prediction law φ be known, with the class prior unknown. Then:
  - (i) the identified set of correctness laws is E(m,p), with m(z) = φ(z)+φ(z̄);
  - (ii) for full-support π, a linear functional is identified iff it has no odd Walsh weight of degree ≥ 3; in particular the accuracy of a positive-weight vote is identified iff the vote is a dictator;
  - (iii) with η known, the identified set is E(m,p) ∩ {TV(π,φ) ≤ 1−η, TV(π, φ∘flip) ≤ η}, and the class-conditional accuracies are identified.
- *Evidence.* Proofs in the earlier notes and Q1. Numerical agreement ≤ 2e−15 (`q1_identification.py`, `q1b_walsh_check.py`, `q5`).
- *Novelty risk:* **medium**. The LP is in Mazzetto 2021; the sign ambiguity is folklore. The Walsh characterisation and the "only dictators" corollary were not found.
- *Feasibility:* **high**. Already proved; 1 week to polish.

**Theorem B (sharp second-order majority-vote bounds) — status PROVED.**
- *Statement.* Given only Gibbs risk μ and tandem risk s (equivalently μ and expected disagreement E[D]):
  - (i) as K → ∞, sup L(MV) = min{C-bound, μ + E[D]}, so L(MV) ≤ E[L] + E[D] always;
  - (ii) for finite odd K and uniform weights, sup L(MV) = min{1, E[ν(r+K−ν)]/(rK), min_j E[(ν−j)(ν−j−1)]/((r−j)(r−j−1))};
  - (iii) for weighted votes, an analogous formula over the value set V_ρ;
  - (iv) with individual (p_i, r_ij), the sharp bound is the pairwise LP. It strictly improves on FO, TND and C, by about 2.2× on real data.
- *Evidence.* `q3_ml_bounds.py`, `q3b_closed_forms.py` (≤ 1.3e−7 vs LP); `q8` real-data table.
- *Novelty risk:* **medium**. (ii) is likely a special case of Boros–Prékopa discrete-moment bounds; the C-bound sharpness region is Germain Prop. 13. "E[L]+E[D]" was not found anywhere.
- *Feasibility:* **high** for the oracle statements. Medium for PAC-Bayes-kl empirical versions (2–3 weeks).

**Theorem C (weights do not help under pairwise ambiguity) — status CONJ, strong NUM evidence.**
- *Statement.* For any (p, r) with p_i > ½, max over positive-weight votes w of min_{π∈P(p,r)} acc_π(w) is attained by the plain majority of some odd subset S ⊆ [K].
  - For K=5 the value is max{p_max, max_{|S|=3} closed form, L(maj₅)}.
  - In the exchangeable case it is 1 − min{e, 3T, 1−5e+10T}.
  - With the full unlabeled law the statement fails: there are explicit K=4 counterexamples with gain 0.067.
- *Evidence.* 0 gain in 2460 adversarial LPs (K=4), 525 adversarial and 376 natural instances (K=5), and 40 random/local-search instances (K=7) (`q6`, `q10`, `q5`). K=5 formula: 300 instances at 1e−15 (`q10`). Symmetric formula: 11,680 cells at 7.5e−15 (`q12`).
- *Novelty risk:* **low–medium**. Not found anywhere. It upgrades Mazzetto's heuristic subset search to max-min-optimal over all weighted votes, and pins down when B&F-style randomisation or weighting can add value.
- *Feasibility:* **medium/uncertain**. K=4 is a finite check over one non-subset type (2,1,1,1) and is doable. General K needs an LP-duality argument on P(p,r) extreme points. Recommendation: prove K ≤ 5 rigorously (computer-assisted via the vertex enumeration of P(p,r)'s dual) and state general K as a conjecture.

**Theorem D (finite-sample certificates and the price of marginal labels) — status PROVED (constants empirical).**
- *Statement.* With p̂ from n_ℓ marginal labels per classifier and m̂ from n_u unlabeled points, with probability ≥ 1−δ, for every self-dual g:
  - L(g) − 2C_K(g)(ε_p + ε_m) ≤ L̂(g) ≤ L(g) ≤ acc(g);
  - ε_p = O(√(log(K/δ)/n_ℓ)) and ε_m = O(√((2^{K−1}+log(1/δ))/n_u)), with ε_m replaced by O(√(K² log(K/δ)/n_u)) for the pairwise set;
  - the certified max-min rule is O(ε_p + ε_m)-optimal;
  - C_K(dictator) = 1, while C_K(maj) is observed to equal K. So a committee is certifiably preferable only if its identified advantage exceeds about (C_K(g)−1)ε_p.
- *Evidence.* Proof in Q3. 100% validity and gap-bound satisfaction in `q7`/`q9`; slope −0.44; real data 120/120 split-ensembles with 100% coverage (`q8`); the crossover at about 3,300 labels per classifier (`q11`).
- *Novelty risk:* **medium**. The technique is standard LP sensitivity, and econometrics has sharper asymptotic inference (Fang et al. 2023). The "price of marginal labels" insight and its link to dictator-selection are new in this ML setting.
- *Feasibility:* **high** for the theorem (1–2 weeks). Proving C_K(maj) ≤ K is medium.

### Inferences
**Six-week plan for 3–4 students.**
- Weeks 1–2: A + D write-up, and a C proof for K=4 and 5.
- Weeks 3–4: B, the PAC-Bayes-kl empirical versions.
- Weeks 5–6: experiments E1–E3 (Q4).

**Biggest risk.** Reviewers may view the work as "Mazzetto 2021 + econometric partial identification". A secondary risk: the empirical story ("the robust rule is usually the best single model") may read as negative. Frame it as **certification** (valid, tight bounds) plus a **structure theorem** (subset majorities suffice; marginal labels force dictators), rather than as a better aggregator.

### Gaps
- The conjecture in C is unproven.
- C_K(g) has no proven closed-form bound.
- No LLM-judge or crowdsourcing experiment was run (blocked hosts).
- Boros–Prékopa equivalence and the novelty of E[L]+E[D] are unverified against full texts.
