# Weighted majority voting when competences are known but dependence is uncertain: what social choice and classifier fusion already know

Scope note. The setting: K binary classifiers/voters, Z_i = 1{i correct}, p_i = P(Z_i = 1) known; the dependence structure is unknown or partly known. "WMR" means weighted majority rule: sign(sum_i w_i s_i), where s_i ∈ {±1} are votes. "Neutral" means the rule treats the two labels symmetrically (self-dual), so it can be written as a function of the correctness pattern Z. Items tagged **[RECENT]** date from 2020 to 2026. Items tagged **[UNVERIFIED]** could not be checked against the source text in this session. Results tagged **[OWN COMPUTATION]** are my own numerical checks, not literature. They used scipy/HiGHS LPs; the scripts are in the session scratchpad, which is not persistent, and the method is described inline so it can be reproduced.

---

## Q1. Optimal weights under independence (Nitzan–Paroush, Shapley–Grofman, extensions)

### Takeaway
Under conditional independence of votes given the truth, and with a uniform prior, the Bayes-optimal rule is a WMR with log-odds weights w_i = log(p_i/(1−p_i)). This is Nitzan & Paroush (1982) and Shapley & Grofman (1984); it is the naive-Bayes rule. With a non-uniform prior, class-asymmetric competences or c > 2 classes, the rule stays a weighted (super)majority with modified weights and threshold. Every one of these results depends on (conditional) independence. None of them says anything about which weights are optimal or robust under dependence.

### Cited Findings
- **Shapley & Grofman (1984)**, "Optimizing group judgmental accuracy in the presence of interdependencies," *Public Choice* 43(3):329–343, DOI 10.1007/BF00118940. Abstract: if individual choices are independent and each alternative has prior 1/2, "the optimal group procedure is a weighted voting rule with weights w_i proportional to log(p_i/(1−p_i))". The paper then "examine[s] the implications for optimal group choice of interdependencies among individual choices". The dependence results themselves could not be read (paywalled) **[UNVERIFIED content]** — [Springer](https://link.springer.com/article/10.1007/BF00118940)
- **Nitzan & Paroush (1982)**, "Optimal decision rules in uncertain dichotomous choice situations," *International Economic Review* 23(2):289–297. Bibliographic data confirmed; no abstract online — [IDEAS](https://ideas.repec.org/a/ier/iecrev/v23y1982i2p289-97.html)
- Exact modern statement of the Nitzan–Paroush theorem (Berend & Kontorovich): assume (i) votes are independent conditional on the truth Y and (ii) P(Y=+1) = P(Y=−1) = 1/2. Then the optimal rule is f_OPT(x) = sign(Σ w_i x_i) with w_i = log(p_i/(1−p_i)), which is naive Bayes. Their error bounds use the "committee potential" Φ = Σ(p_i − 1/2)w_i: P(err) ≤ exp(−Φ/2) and P(err) ≥ 3/(8[1+exp(2Φ+4√Φ)]). They call independence "the truly restrictive assumption", and the paper does not analyse dependence — Berend & Kontorovich, "Consistency of weighted majority votes," NeurIPS 2014, [PDF](https://proceedings.neurips.cc/paper_files/paper/2014/file/3cf61323ed04cc4581f435dab1256ab8-Paper.pdf); [arXiv 1312.0451](https://arxiv.org/abs/1312.0451)
- Open problems listed by Berend & Kontorovich (2014), both for *estimated* competences under independence: a better sample-size dependence for the adaptive bound, and an error estimate for the Bayesian empirical rule — [PDF](https://proceedings.neurips.cc/paper_files/paper/2014/file/3cf61323ed04cc4581f435dab1256ab8-Paper.pdf)
- **Ben-Yashar & Nitzan (1997)**, "The optimal decision rule for fixed-size committees in dichotomous choice situations: the general result," *International Economic Review* 38(1):175–186. Abstract: "derives the optimal decision rule for fixed-size committees" and "generalizes a number of earlier results" — [IDEAS](https://ideas.repec.org/a/ier/iecrev/v38y1997i1p175-86.html). **[UNVERIFIED exact statement]** My recollection, not checked against the text: with state-dependent competences α_i = P(correct | state 1) and β_i = P(correct | state 2) and prior π, the optimal rule is the likelihood-ratio (naive-Bayes) rule. That is a weighted rule with weights log[α_i β_i / ((1−α_i)(1−β_i))]/2 and a threshold that depends on log(π/(1−π)) and Σ log(α_i/β_i) asymmetries, i.e. a weighted *qualified* majority.
- **Multi-class version** (Kuncheva & Rodríguez 2014, *Knowledge and Information Systems* 38:259–275, DOI 10.1007/s10115-012-0586-6). Assumptions: "the classifiers give their decisions independently conditioned upon the class label", accuracy p_i the same for every class, and errors spread uniformly over the other c−1 classes. Under these, log P(ω_k|s) ∝ log P(ω_k) + Σ_{i votes k} [log(p_i/(1−p_i)) + log(c−1)], so the per-vote weight is log((c−1)p_i/(1−p_i)). The paper builds a hierarchy MV ⊂ WMV ⊂ Recall combiner ⊂ Naive Bayes by "progressively relaxing" the accuracy assumption. Independence is never relaxed: "it is not clear how the violation of the assumption of conditional independence affects the performance." It also notes that "the optimality of the combiner … is asymptotic" (it holds as sample size → ∞) — [PDF](https://lucykuncheva.co.uk/papers/lkjrKAIS14.pdf)
- **Baharad & Ben-Yashar (2009)**, "The robustness of the optimal weighted majority rule to probability distortion," *Public Choice* 139:53–59. Under Tversky–Kahneman probability distortion, "neither the rule nor the weights are changed". This is "robustness" in a behavioural sense, not with respect to dependence — [Springer](https://link.springer.com/doi/10.1007/s11127-008-9378-7)
- **Ben-Yashar & Danziger (2015)**, "When is voting optimal?" *Economic Theory Bulletin* 3(2):341–356. The optimal rule depends only on members' posteriors. The conditions under which it coincides with weighted (super)majority voting are "very stringent", and "voting is generally not an efficient way to make collective choices." Whether dependence is allowed is not visible from the abstract page — [Springer](https://link.springer.com/article/10.1007/s40505-014-0065-1)
- **Berend & Sapir (2001)**, "Optimality of the expert rule under partial information," *Acta Applicandae Mathematicae* 69:141–162. Competences are i.i.d. random draws, only the ranking is known, and choices are independent. It gives an explicit formula for P_e(n), the probability that the expert rule is the optimal (log-odds) rule. Example: P_e(3) = 3·(1/2)² = 0.75 for exponentially distributed log-expertise — [PDF](https://tzin.bgu.ac.il/~lsapir/Papers/erule.pdf)
- **Berend & Sapir (2003)**, "Between the expert and majority rules," *Advances in Applied Probability* 35(4):941–960. Computes optimality probabilities for "restricted majority rules" and "balanced expert rules" (families between the expert and majority rules) — [Cambridge](https://www.cambridge.org/core/journals/advances-in-applied-probability/article/between-the-expert-and-majority-rules/E1DFDAB8995FD3A65B3B6A0650069430)
- **[RECENT]** Romaniega Sancho, "Unweighted Condorcet jury theorem and miracle of aggregation do not hold almost surely" (arXiv 2108.00733). Under a centred prior on competences, the CJT fails almost surely. Weighted majority with competence-correlated weights (each weight ≥ 1) restores it almost surely. Independence is assumed; correlated votes appear only as future work — [arXiv](https://arxiv.org/pdf/2108.00733v3)

### Inferences
- The log-odds rule is simply the Bayes posterior when the likelihood factorises. Any departure from conditional independence changes the likelihood, so the Bayes-optimal rule generally stops being a WMR. That is why the dependent case cannot be handled by re-weighting alone (see Q2, Jaffe et al.).
- In the project's notation, the independence-optimal weights are w_i = log(p_i/(1−p_i)). They are the natural baseline to compare against any robust weights.

### Gaps
- Exact theorem statements of Nitzan & Paroush (1982) and Ben-Yashar & Nitzan (1997) could not be read (JSTOR/paywall). The statement above comes from Berend & Kontorovich's restatement.
- Shapley & Grofman's results under interdependence are not visible; the abstract only says they are "examined".

---

## Q2. Optimal weights or decision rules when votes are correlated with a known structure

### Takeaway
The correlated-jury literature is overwhelmingly about **simple majority** (or a quota), studying how correlation changes its accuracy. I found no paper that derives optimal *weights* as a function of competences and a given correlation matrix for heterogeneous voters. The closest results are:
- Kaniovski & Zaigraev's correlation-robust *quota* for homogeneous juries.
- Berg (1994), which evaluates (does not optimise) some WMRs under dependence.
- Nitzan & Paroush (1984), where the optimal rule is identified only for "total interdependence".
- Jaffe et al. (2016), where the optimal rule under latent-group dependence is non-linear, so it is not a WMR.

### Cited Findings
**Homogeneous juries (exchangeable votes): Kaniovski–Zaigraev programme**
- **Kaniovski & Zaigraev (2011)**, "Optimal jury design for homogeneous juries with correlated votes," *Theory and Decision* 71(4):439–459 ([IDEAS](https://ideas.repec.org/a/kap/theord/v71y2011i4p439-459.html)). Read from the 2009 preprint ([PDF](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_robust.pdf)); theorem numbering may differ in the published version. Common competence p and common pairwise correlation c, with a Bahadur distribution truncated at second order. Results:
  - Theorem 1: the probability of at least k correct votes is P^k_n(p,c) = P^k_n(p,0) + 0.5·c·n(n−1)·((k−1)/(n−1) − p)·C(n−1,k−1)·p^{k−1}(1−p)^{n−k}.
  - Corollary 1: P^k_n decreases in c for k < p(n−1)+1 and increases for k > p(n−1)+1.
  - The **"correlation-robust voting quota"** α = p(n−1)+1 makes accuracy independent of c. When α is not an integer, "⌊α⌋ [is] the optimal correlation-robust voting rule." The authors concede that it "would not necessarily maximize" accuracy.
  - Theorem 2, n = 3: M_3(p,c) − p = p(1−p)(2p−1)(1−3c). Majority beats a single juror iff c < 1/3.
  - Unanimity is always worse than a single juror.
  - The paper has **no minimax criterion, no LP, and no weights.**
- **Zaigraev & Kaniovski (2012)**, "Bounds on the competence of a homogeneous jury," *Theory and Decision* 72(1):89–112 ([IDEAS](https://ideas.repec.org:443/a/kap/theord/v72y2012i1p89-112.html)); preprint ([PDF](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_minmax.pdf)). Exchangeable votes with all correlations of every order free. Over that set, simple majority (odd n) has min M = (n(2p−1)+1)/(n+1) and max M = min{2np/(n+1), 1}, computed by LP. For n = 3: [(3p−1)/2, min{1.5p, 1}]; at p = 0.75 this is [0.625, 1] against 0.844 under independence. Simple majority only.
- **Kaniovski (2009)**, "An invariance result for homogeneous juries with correlated votes," *Mathematical Social Sciences* 57(2):213–222. Defines "second-order invariant" joint distributions, under which Condorcet's probability does not depend on pairwise correlations, and shows they exist for any jury size — [IDEAS](https://ideas.repec.org/a/eee/matsoc/v57y2009i2p213-222.html)

**Heterogeneous competences with given pairwise correlations**
- **Kaniovski (2010)**, "Aggregation of correlated votes and Condorcet's Jury Theorem," *Theory and Decision* 69(3):453–468, DOI 10.1007/s11238-008-9120-4 ([IDEAS](https://ideas.repec.org:443/a/kap/theord/v69y2010i3p453-468.html)); preprint ([PDF](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet.pdf)).
  - Given p_i and c_ij, it computes LP min/max of simple-majority accuracy.
  - It also defines a QP "Q-solution", the feasible distribution closest to independence, whose majority accuracy equals the independent value (Theorem 2), and an entropy "E-solution".
  - Table 2, n = 3:
    - p = 0.75, c = 0: bounds [0.813, 0.938].
    - p = 0.75, c = 0.2: bounds [0.70, 0.90].
    - Heterogeneous p = (0.75, 0.6, 0.6), c = 0.2: bounds [0.577, 0.792].
  - **No weighted voting.**
- Kaniovski's 2010 slides summarise the programme and state no results on weighted voting or optimal weights. He states that "Empirical evidence overwhelmingly refutes the assumption of independent votes" (US Supreme Court data) — [LSE slides](https://www.lse.ac.uk/cpnss/assets/documents/voting-power-and-procedures/workshops/2011/Kaniovski-Slides.pdf)

**Weighted rules under dependence (classical)**
- **Nitzan & Paroush (1984)**, "The significance of independent decisions in uncertain dichotomous choice situations," *Theory and Decision* 17(1):47–60, DOI 10.1007/BF00140055. Abstract: "independent voting is always weakly superior to any pattern of interdependence among individual decisions"; for "total interdependence" "the optimal rule is explicitly identified" — [Springer](https://link.springer.com/article/10.1007/BF00140055). **[UNVERIFIED scope]** "Pattern of interdependence" probably means deterministic following or copying structures. It cannot mean arbitrary correlation, because Kaniovski (2010, Cor. 1.1) shows that negative correlation *raises* the accuracy of a homogeneous competent jury ([preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet.pdf)). Check the full text before citing.
- **Berg (1994)**, "Evaluation of some weighted majority decision rules under dependent voting," *Mathematical Social Sciences* 28(2):71–83. No abstract online — [IDEAS](https://ideas.repec.org/a/eee/matsoc/v28y1994i2p71-83.html). Karge et al. (2024) describe it as extending the Boland et al. dependence models to dichotomous weighted voting rules ([AAMAS'24 PDF](https://iccl.inf.tu-dresden.de/w/images/e/e3/AAMAS24_OL.pdf)). **This is the single closest classical paper to "WMRs under dependence"; it evaluates rules rather than optimising them, and its content is [UNVERIFIED].**
- **Boland, Proschan & Tong (1989)**, "Modelling dependence in simple and indirect majority systems," *Journal of Applied Probability* 26(1):81–88. Introduces two dependence models (including the "opinion leader" model), proves monotonicity results, and compares direct with indirect majority systems — [Cambridge](https://www.cambridge.org/core/journals/journal-of-applied-probability/article/abs/modelling-dependence-in-simple-and-indirect-majority-systems/070D6335BDDDDC7AF4D70BC9B21B0B7B)
- Other correlated-CJT references, as listed in Zaigraev & Kaniovski's bibliography ([preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_minmax.pdf)):
  - Ladha (1992) *AJPS* 36:617–634
  - Ladha (1993) *Social Choice and Welfare* 10:69–85, de Finetti/exchangeability
  - Ladha (1995) *JEBO* 26(3):353–372 ([IDEAS](https://ideas.repec.org/a/eee/jeborg/v26y1995i3p353-372.html))
  - Berg (1993) *SCW* 10:87–95
  - Berg (1985) *Public Choice* 47:377–387
  - Boland (1989) *The Statistician* 38:181–189
  - Berend & Sapir (2007) *SCW* 28:507–528, monotonicity with dependent voters

  All concern simple majority, to the extent visible.

**General dependent juries (asymptotics, simple majority)**
- **Peleg & Zamir (2012)**, "Extending the Condorcet Jury Theorem to a general dependent jury," *Social Choice and Welfare* 39(1):91–125 ([IDEAS](https://ideas.repec.org/a/spr/sochwe/v39y2012i1p91-125.html)); DP 477 ([PDF](https://ratio.huji.ac.il/files/dp477.pdf)).
  - Sufficient condition: (p̄_n − 1/2)/√E(X̄_n − p̄_n)² → ∞. In particular, non-positive pairwise covariances plus √n(p̄_n − 1/2) → ∞ suffice.
  - Necessary-type conditions; non-CJT distributions form a convex set that is dense.
  - Simple majority only ("determined by a simple majority"); no weighted rules; no finite-n results.
- **Pivato (2017)**, "Epistemic democracy with correlated voters," *Journal of Mathematical Economics* 72:51–69. A "suitably chosen voting rule" converges to the truth if average correlation vanishes; covers social-network and DeGroot models. Its bibliography includes Berg (1994) and Pivato (2013) "Voting rules as statistical estimators". Whether it optimises weights is **[UNVERIFIED]** — [IDEAS](https://ideas.repec.org/p/pra/mprapa/69546.html)
- **Dietrich & Spiekermann**, "Independent opinions? On the causal foundations of belief formation and jury theorems" (*Mind*, 2013; vol/pages [UNVERIFIED]). Argues that the right independence notion is conditional on the problem and its common causes — [PhilArchive](https://philarchive.org/rec/DIEIO). Not fetched in detail.

**Opinion-leader models and recent dependence work**
- **[RECENT] Karge, Burkhardt, Rudolph & Rusovac (AAMAS 2024)**, "To lead or to be led: a generalized CJT under dependence." Opinion leader with competence p̂; each agent follows the leader with probability π. Asymptotic CJT for approval voting holds iff π < Δp/(Δp+1); finite-n sample-size bounds. **No weighted rules, no small-n results, no comparison with the best individual** — [PDF](https://iccl.inf.tu-dresden.de/w/images/e/e3/AAMAS24_OL.pdf)
- **[RECENT] Karge (arXiv 2609.00010, July 2026)**, "Bridging dependence: from opinion leader influence to standard correlation measures." Translates the opinion-leader model into covariance and correlation-coefficient terms to unify dependence models — [arXiv](https://arxiv.org/abs/2609.00010)
- **[RECENT] Mulligan (2024)**, "Optimizing political influence: a jury theorem with dynamic competence and dependence," *Social Choice and Welfare* 63(3):509–530. Simple majority. Notes that competence plus pairwise correlation "do not determine a unique joint distribution" — [Springer](https://link.springer.com/article/10.1007/s00355-022-01407-5)
- **[RECENT] Braha & de Aguiar (arXiv 2607.14288, July 2026)**, majority correctness after network interaction. Contrarian updating drives accuracy to 1/2; topology shapes vote correlations — [arXiv](https://arxiv.org/abs/2607.14288)
- Thiele (2017), "A note on the Condorcet jury theorem for couples," *Theory and Decision* 83:355–364. Couple-wise dependence, simple majority — [Springer](https://link.springer.com/article/10.1007/s11238-017-9602-3)

**ML side: known dependence structure**
- **Jaffe, Fetaya, Nadler, Jiang & Kluger (AISTATS 2016, PMLR 51)**, "Unsupervised ensemble learning with dependent classifiers." Latent-group model: classifiers in a group share a latent variable, and groups are conditionally independent given Y. "The optimal ensemble learner is in general non-linear in the m classifiers"; the linear (log-odds) rule "may be suboptimal, even when provided with the true classifier accuracies." Finding the groups from a sample covariance is NP-hard in general (Lemma 3) — [arXiv](https://arxiv.org/pdf/1510.05830); [PMLR](https://proceedings.mlr.press/v51/jaffe16.html)

### Inferences
- With full knowledge of the joint distribution π over {0,1}^K, the optimal *neutral* rule is trivial: for each complementary pair (z, 1−z), output the member with the larger π. Such a rule is generally not a WMR, and for K ≥ 6 it may not even be monotone or threshold. This is why the literature studies quotas, invariance, or bounds rather than "optimal weights under correlation".
- Kaniovski & Zaigraev's correlation-robust quota is the only "robust rule choice under correlation" result in the CJT literature I found. It is invariance-based (accuracy does not depend on c), not maxmin, and is restricted to the second-order Bahadur family. The truncated Bahadur family is *not* the worst case, as the gap between Table 1 and the LP bounds in Kaniovski 2010 shows.

### Gaps
- No source found that gives optimal WMR weights as a closed-form function of (p_i, ρ_ij) for heterogeneous voters, even for K = 3. Berg (1994) is the most likely candidate and should be read in full.
- Romeijn & Atkinson and Kanazawa: not located in this session. A 1998 *MSS* note "A brief note on a further refinement of the CJT for heterogeneous groups" (35(1):69–73) appeared in search results ([IDEAS](https://ideas.repec.org/a/eee/matsoc/v35y1998i1p69-73.html)); attributing it to Kanazawa is [UNVERIFIED], and it appears to be about independent voters.

---

## Q3. Robust (maxmin, worst-case or ambiguity) weight choice when dependence is unknown or partly known

### Takeaway
1. **Marginals only.** The question is *solved*, and the answer is negative for weighting. In binary-state, binary-action problems, when only each source's marginal informativeness is known and the correlation is adversarial, following the single best source is robustly (maxmin) optimal, and no aggregation rule, weighted or otherwise, beats it. This is de Oliveira, Ishii & Lin, Theorem 1. With symmetric agents, a random dictator is minimax- and regret-optimal (Arieli et al.).
2. **Partial dependence information.** The non-trivial robust problem needs some dependence information: pairwise joint accuracies or agreements, correlation bounds, or the full distribution of vote profiles.
   - The ML paper closest to the project is Balsubramani & Freund (COLT 2015). They solve the minimax game when the *joint distribution of predictions* (from unlabelled data) is known and only lower bounds on accuracies are known. The minimax-optimal predictor is a clipped weighted vote with weights from a convex program.
   - In economics, Levy & Razin (2022) handle forecasts with a bounded-correlation ambiguity set.
   - I found **no paper** that chooses WMR weights by maxmin over joint distributions consistent with known marginals plus known pairwise joint success probabilities (or pairwise correlations) for K = 5, or general K. In the social-choice/CJT literature I found no maxmin weight choice at all.

### Cited Findings
**Marginals only: best source / random dictator is robust-optimal**
- **[RECENT] de Oliveira, Ishii & Lin**, "Robust Aggregation of Correlated Information" (Sept. 2026 version); earlier "Robust Merging of Information" (ACM EC 2021, DOI 10.1145/3465456.3467533; arXiv 2106.00088).
  - Model: the decision maker knows each marginal experiment P_j, does not know the joint distribution, and maximises the worst-case payoff over all joint experiments consistent with the marginals.
  - **Theorem 1:** "For all (A, u) with |A| = |Θ| = 2, any best-source strategy is robustly optimal", and V(P_1,…,P_m) = max_j V(P_j). Uniqueness: "if each marginal experiment has full support" and the best source is unique, "then all robustly optimal strategies are best-source strategies."
  - Further results: with more actions, decompose into binary comparisons and use at most |A|−1 sources (Theorem 2, Corollary 2). An added source helps only if it is best in some binary subproblem (Corollary 3). With general states, the value is a concavification (Theorem 3).
  - Sources: [PDF](https://papers.xiao-lin.com/Robust_Information_Aggregation.pdf); [arXiv](https://arxiv.org/pdf/2106.00088); [ACM](https://dl.acm.org/doi/10.1145/3465456.3467533)
- **[RECENT] Arieli, Babichenko, Talgam-Cohen & Zabarnyi**, "Universally robust information aggregation for binary decisions," also circulated as "A random dictator is all you need" (arXiv 2302.03667). Symmetric agents with binary recommendations; the decision maker knows the prior and each agent's marginal, and the correlation is adversarial.
  - For many agents, random dictator is the unique optimal anonymous rule under minimax, regret and approximation-ratio criteria. Minimal regret for any n is characterised by concavification.
  - Simple majority is "far from optimal".
  - In the non-symmetric case, de Oliveira et al. show "follow the best agent" is minimax-optimal.
  - Sources: [ar5iv](https://ar5iv.labs.arxiv.org/html/2302.03667); [arXiv v2](https://arxiv.org/abs/2302.03667v2). Venue **[UNVERIFIED]**; a search listing shows DOI 10.1145/3580507.3597710, i.e. ACM EC 2023.
- **[OWN COMPUTATION]** Consistency check in the project's own setting, where only P(Z_i = 1) = p_i is known, less than the full marginal experiment. For 200 random p ∈ [0.5, 0.99]^K, I computed the LP worst-case accuracy of every WMR (all 4 labelled WMRs for K = 3 and all 81 for K = 5) over all joint distributions with those marginals. The best worst case among WMRs equalled max_i p_i to within 2·10⁻¹⁶ in every instance, so dictatorship of the best classifier is robust-optimal. This agrees with de Oliveira et al. and with the homogeneous min M = (n(2p−1)+1)/(n+1) < p of Zaigraev & Kaniovski ([preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_minmax.pdf)).

**Partial dependence information**
- **Balsubramani & Freund (2015)**, "Optimally combining classifiers using unlabeled data," COLT 2015, PMLR 40:211–225.
  - Setup: the prediction matrix F of the ensemble on n unlabelled points is known, and so are lower bounds b on each classifier's label correlation, (1/n)Fz ≥ b. The adversary picks labels z ∈ [−1,1]^n.
  - Theorem 2: the value is V = −γ(σ*), with slack function γ(σ) = (1/n)Σ_j[|x_jᵀσ| − 1]₊ − bᵀσ minimised over σ ≥ 0 (convex). The optimal prediction is g_j* = clip(x_jᵀσ*, [−1,1]), a clipped weighted vote that randomises inside the "zero box".
  - If F is unknown, the solution reduces to "the classifier with the best error guarantee". With F known, a weighted combination "can perform significantly better than any single classifier"; in an extremal example the gain is p-fold.
  - Sources: [PMLR](https://proceedings.mlr.press/v40/Balsubramani15.html); [arXiv](https://arxiv.org/pdf/1503.01811)
- **Arachie & Huang (2021)**, "A general framework for adversarial label learning," *JMLR* 22:1–33. Weak signals with expected-error bounds b_i; the adversary picks labels within the bounds; the learner minimises the worst-case error (the inner problem is an LP). Correlated *copies* of a signal leave the optimum unchanged, unlike majority vote, where redundant signals "outvote". Beating the best weak signal is shown only empirically — [JMLR PDF](https://jmlr2020.csail.mit.edu/papers/volume22/20-537/20-537.pdf)
- **Mazzetto, Sam, Park, Upfal & Bach (AISTATS 2021, PMLR 130:3196–3204)**: no independence assumption; majority votes over selected small subsets of sources, with provable error — [PMLR](https://proceedings.mlr.press/v130/mazzetto21a.html). **Mazzetto, Cousins, Sam, Bach & Upfal (ICML 2021, PMLR 139:7534–7543)**: adversarial labelings constrained by estimated source statistics, sources "may be arbitrarily correlated" — [PMLR](https://proceedings.mlr.press/v139/mazzetto21a.html)
- **Levy & Razin (2022)**, "Combining forecasts in the presence of ambiguity over correlation structures," *Journal of Economic Theory* 199:105075, DOI 10.1016/j.jet.2020.105075. Decision makers "understand each information source in isolation but are uncertain about the correlation". The set of predictions is pinned down by the naive-Bayes (correlation-neglect) interpretation and a bound on correlation measured by pointwise mutual information. A larger bound gives "more conservative behaviour". This is about forecasts and posteriors, not voting weights; the exact decision criterion was not visible — [IDEAS](https://ideas.repec.org/a/eee/jetheo/v199y2022ics0022053118303144.html)
- **[RECENT] Pan, Chen & Kong**, "Robust decision aggregation with second-order information" (WWW 2024, DOI 10.1145/3589334.3645384; arXiv 2311.14094). **Two experts only.** The aggregator sees the recommendations plus each expert's prediction of the other's recommendation, and is scored by minimax regret.
  - With general information structures, regret is ≥ 0.5 for every aggregator, and second-order information "provides no benefit".
  - With conditionally independent signals: a threshold aggregator reaches 1/3 (deterministic), the uniform aggregator reaches 0.25, and in the homogeneous case the value is 3 − 2√2.
  - Open: more experts; closing the 0.1667 vs 0.1673 gap.
  - Sources: [arXiv](https://arxiv.org/html/2311.14094v1)
- **[RECENT] Kong**, "Peer expectation in robust forecast aggregation: the possibility/impossibility" (arXiv 2402.06062). Two forecasters, conditionally independent or refinement-ordered — [arXiv](https://arxiv.org/html/2402.06062v1)

**"Robust voting" papers that are *not* about competence and dependence (false positives to avoid)**
- **[RECENT] Nakada, Nitzan & Ui (arXiv 2507.22655, July 2025)**, "Robust voting under uncertainty." Uncertainty is over the joint distribution of *preference* profiles, which may be arbitrarily correlated. A rule is P-robust iff ∃ w ≥ 0 such that the weighted-average responsiveness exceeds 1/2 for every p ∈ P. When P is all distributions, the robust rules are exactly **WMRs without ties**; for odd n, simple majority is the unique anonymous robust rule. The criterion is Pareto-undominance, not maxmin. There are no competences, but it gives a characterisation of WMRs via robustness that may be useful as a tool — [arXiv](https://arxiv.org/html/2507.22655)
- **[RECENT] Kiguchi, Noguchi & Takahashi (2026)**, "Robustly optimal voting rule," *Journal of Economic Theory* 236:106219. A utilitarian, strategy-proof mechanism with worst-case welfare or regret; the median voter rule is uniquely optimal. Not epistemic — [summary](https://eurekamag.com/research/107/915/107915857.php); [author page](https://sites.google.com/view/noriaki-kiguchi/home/)

**[OWN COMPUTATION] Structure of the K = 3 robust problem with known pairwise information**
- Derivation, elementary: for K = 3, fixing (p_1, p_2, p_3, p_12, p_13, p_23) and normalisation gives 7 independent linear constraints on 8 cell probabilities. The consistent joint distributions therefore form a segment parametrised by t = P(Z = 111):
  - π_t = π_0 + t·d with d(z) = (−1)^{3−|z|};
  - t ∈ [max(0, p_12+p_13−p_1, p_12+p_23−p_2, p_13+p_23−p_3), min(p_12, p_13, p_23, 1 − Σp_i + Σp_ij)].
- Every neutral rule h has accuracy A(h) + t·Δ(h) with Δ(h) = Σ_z h(z)d(z) ∈ {0, ±2, ±4}:
  - simple majority: Δ = −2;
  - the 3 dictators: Δ = 0, i.e. invariant to third-order dependence;
  - "signed" votes such as sign(s1 + s2 − s3): Δ = +2.
- Check: homogeneous p = 0.75 with pairwise independence (p_ij = 0.5625) gives t ∈ [0.375, 0.4375] and majority accuracy 1.6875 − 2t ∈ [0.8125, 0.9375]. This reproduces Kaniovski 2010 Table 2 ([0.813, 0.938], [preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet.pdf)). Here majority's worst case (0.8125) beats the best single classifier (0.75).
- Numerical experiment: random true π ~ Dirichlet(0.7·1) on {0,1}³, 184 instances with all p_i ≥ 0.5, all 16 deterministic neutral rules evaluated.
  - In **22 of 184 instances**, the best rule by worst case was a **signed WMR with a negative weight** on a classifier whose p_i ≥ 0.5. The winning weight vectors were (1,1,−1), (1,−3,3) and (−3,1,3). They beat every monotone WMR by up to **0.112**.
  - **Randomised** neutral rules (maxmin LP via duality) never beat the best deterministic neutral rule; the gap was ≤ 2·10⁻¹⁶. This is a numerical observation, not a proof.
  - In a separate run of 2000 Dirichlet instances without the p ≥ 0.5 filter, majority's worst case exceeded max_i p_i in 3.0% of instances, by up to 0.152. This frequency depends on the sampling distribution.

### Inferences
- **Equivalence (my inference, worth verifying carefully).** For K = 3, knowing the pairwise joint success probabilities p_ij given the p_i is equivalent to knowing the pairwise *agreement* rates, since P(Z_i = Z_j) = 1 − p_i − p_j + 2p_ij. That in turn is equivalent to knowing the full distribution of vote profiles up to a global flip: the 4 complementary-pair masses are an invertible function of the 3 agreement rates. The K = 3 "marginals + pairwise" robust problem is therefore the population, label-symmetric version of Balsubramani & Freund's game, with **equality** accuracy constraints in place of lower bounds. Two consequences:
  - (a) Agreement rates can be estimated from unlabelled data, which gives a good ML motivation.
  - (b) In Balsubramani & Freund, lower-bound constraints force σ ≥ 0 (nonnegative weights). With exact p_i, the dual variables are sign-free, which explains the negative-weight robust rules found above.
- For K ≥ 4 the equivalence breaks. Pairwise information (K(K−1)/2 numbers) no longer determines the distribution of vote profiles (2^{K−1} − 1 numbers; 15 for K = 5 against 10 pairwise). The project's "pairwise-only" game is then strictly harder for the decision maker than Balsubramani & Freund's game, because the adversary also chooses the higher-order agreement structure. I found no paper analysing this game.
- Kaniovski & Zaigraev's correlation-robust quota and the K = 3 rules with Δ(h) = 0 are both instances of "dependence-invariant rules". Generalising this (rules whose accuracy is determined by moments up to order m) gives a clean, testable notion.

### Gaps
- No found work gives **maxmin-optimal WMR weights given marginals plus pairwise information for K = 5** (or general K), nor characterises when the robust-optimal rule is a (signed) WMR or how its weights relate to log-odds.
- No found work studies ambiguity sets of the form "|ρ_ij − ρ̂_ij| ≤ ε" (or bounds on correlations) for *voting* rules; Levy & Razin do this for forecasts with a PMI bound.
- Exact venue of Arieli et al. and the full theorem statements (Theorems 4.3 and 5.1) were not verified.

---

## Q4. Classifier-fusion literature: weighted majority vote, optimal weights, bounds, beating the best classifier

### Takeaway
Classifier fusion inherits the Nitzan–Paroush/Shapley–Grofman log-odds weights (Kuncheva's book §4.4.2; Kuncheva & Rodríguez 2014), always under conditional independence. Dependence is handled only through *bounds* on simple-majority accuracy:
- Matan; Kuncheva et al. 2003 (patterns of success and failure);
- Ruta & Gabrys 2002;
- Narasimhamurthy's LP bounds;
- PAC-Bayesian C-bound and tandem bounds, which involve pairwise disagreement for a fixed weighting ρ.

The PAC-Bayes line (MinCq; tandem-bound minimisation) *does* choose weights using second-order (pairwise) information, but it minimises an upper bound, not the exact worst case.

### Cited Findings
- **Kuncheva, *Combining Pattern Classifiers*, 2nd ed., Wiley 2014.** From the table of contents (body not accessible):
  - §4.3.3 "Limits on the Majority Vote Accuracy: An Example";
  - §4.3.4 "Patterns of Success and Failure";
  - §4.3.5 "Optimality of the Majority Vote Combiner";
  - §4.4 "Weighted Majority Vote", with §4.4.2 "Optimality of the Weighted Majority Vote Combiner";
  - Appendix 4.A.1 "Matan's Proof for the Limits on the Majority Vote Accuracy".

  The preface says Ch. 4 defines optimality conditions "by progressively relaxing an assumption" — [PDF (front matter)](https://lucykuncheva.co.uk/Combining_Pattern_Classifiers_Methods_and_Algorithms_2nd_ed_Kuncheva%202014-09-09.pdf). Course slides based on the book give b_i ∝ log(p_i/(1−p_i)) and warn "Minimum error not guaranteed" — [RIT slides](https://www.cs.rit.edu/~rlaz/prec2010/slides/MajorityVote.pdf)
- **Kuncheva & Rodríguez (2014)**, *KAIS* 38:259–275. Exact assumptions and formulas are given under Q1. Experiments on 73 datasets found no clear winner among MV, WMV, Recall and NB; NB was slightly better, and WMV did better with few, unbalanced classes. Under accurate estimates *and* independence the ranking is NB > REC > WMV > MV — [PDF](https://lucykuncheva.co.uk/papers/lkjrKAIS14.pdf)
- **Kuncheva, Whitaker, Shipp & Duin (2003)**, "Limits on the majority vote accuracy in classifier fusion," *Pattern Analysis and Applications* 6(1):22–31 — [Springer](https://link.springer.com/doi/10.1007/s10044-002-0173-7). Already covered in the project; patterns of success and failure.
- **Ruta & Gabrys (2002)**, "A theoretical analysis of the limits of majority voting errors for multiple classifier systems," *Pattern Analysis and Applications* 5(4):333–350. Negatively correlated classifiers improve MV, but "negatively correlated classifiers represent a very unlikely situation in real-world classification problems". Multistage organisations widen the limits — [Bournemouth](https://staffprofiles.bournemouth.ac.uk/display/journal-article/5054)
- **Narasimhamurthy (2005)**, "Theoretical bounds of majority voting performance for a binary classification problem," *IEEE TPAMI* 27(12):1988–1995, and the 2003 LNCS 2749 framework paper. LP bounds for MV with pairwise information. Citation taken from the project bibliography; **not re-verified in this session**.
- **Lam & Suen**: cited in Kuncheva & Rodríguez as the 1997 analysis of odd versus even L and of adding or removing classifiers, plus "Optimal combination of pattern classifiers" (1995) — [PDF](https://lucykuncheva.co.uk/papers/lkjrKAIS14.pdf). Full citation **[UNVERIFIED]**, from memory: Lam & Suen (1997), *IEEE Trans. SMC-A* 27(5):553–568; Lam & Suen (1995), *Pattern Recognition Letters* 16:945–954.
- **Germain, Lacasse, Laviolette, Marchand & Roy (2015)**, "Risk bounds for the majority vote: from a PAC-Bayesian analysis to a learning algorithm," *JMLR* 16:787–860.
  - **C-bound**: R(B_Q) ≤ 1 − (1 − 2R(G_Q))²/(1 − 2d_Q), where R(G_Q) is the Gibbs risk and d_Q is the expected pairwise disagreement. It follows from Cantelli–Chebyshev applied to the margin.
  - Tightness: for any margin distribution with 0 < μ₂ ≤ μ₁, a two-point margin distribution attains it (Prop. 13). It is optimal iff C_Q ≤ 2R(G_Q) (Prop. 14).
  - **MinCq** chooses weights by minimising the C-bound, which "reduces to a simple quadratic program".
  - Sources: [ar5iv](https://ar5iv.labs.arxiv.org/html/1503.08329)
- **Masegosa, Lorenzen, Igel & Seldin (NeurIPS 2020)**, "Second order PAC-Bayesian bounds for the weighted majority vote."
  - Tandem bound L(MV_ρ) ≤ 4·E_{ρ²}[L(h,h′)], where L(h,h′) is the joint-error probability. First-order bound L(MV_ρ) ≤ 2·E_ρ[L(h)].
  - No independence assumption. The C-bound is always at least as tight as the tandem bound but harder to estimate and optimise.
  - Minimising the tandem bound to choose weights does not degrade test error; minimising the first-order bound does.
  - Sources: [arXiv](https://arxiv.org/pdf/2007.13532); [NeurIPS](https://papers.nips.cc/paper/2020/hash/386854131f58a556343e056f03626e00-Abstract.html)
- **[RECENT] Cobbenhagen, Carè, Campi, Ramponi, Antunes & Heemels (2021)**, "Novel bounds on the probability of misclassification in majority voting: leveraging the majority size," *IEEE Control Systems Letters* 5(5). Weighted voters, no independence assumption. "The factor of 2 in the bound [B_Q ≤ 2G_Q] was shown to be tight", so MV can be worse than the average voter. Gives conditional bounds given the agreement level — [PDF](https://marco-campi.unibs.it/pdf-pszip/majority%20voting.pdf)
- **Fumera & Roli (2005)**, "A theoretical and experimental analysis of linear combiners for multiple classifier systems," *IEEE TPAMI* (2005). Weighted *averaging* of soft outputs with correlated estimation errors. Listing only; content and volume/pages **[UNVERIFIED]** — [IRIS UniCa](https://iris.unica.it/handle/11584/97944)
- **Can a weighted vote beat the best single classifier in the worst case over dependence?**
  - With marginals only: **no**. See Q3: de Oliveira et al. Theorem 1 ([PDF](https://papers.xiao-lin.com/Robust_Information_Aggregation.pdf)); Balsubramani & Freund's reduction to the best single classifier when F is unknown ([arXiv](https://arxiv.org/pdf/1503.01811)); and the own LP check over all K = 3 and K = 5 WMRs.
  - With the joint prediction pattern known: **yes, possibly by a large factor** (Balsubramani & Freund).
  - With pairwise information at K = 3: yes in some instances; e.g. Kaniovski 2010's p = 0.75, c = 0 case has majority worst case 0.813 > 0.75 ([preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet.pdf)).

### Inferences
- The fusion literature's "optimal weights" are entirely the independence result. Its dependence results are bounds for *simple* MV (Matan/Kuncheva/Narasimhamurthy) or *upper bounds* for a fixed ρ (C-bound/tandem). Exact worst-case weight optimisation under pairwise constraints sits in the gap between the two.
- MinCq and tandem-bound minimisation are natural ML baselines for any "robust weights from pairwise information" paper. The C-bound is a relaxation, while the LP worst case is exact; a comparison would show how loose the C-bound weights are for small K.

### Gaps
- Kuncheva's §4.4.2 statement (exact assumptions, and whether the c-class log(c−1) term is included) could not be read; Kuncheva & Rodríguez 2014 is the verified substitute.
- Narasimhamurthy 2005 and Lam & Suen 1997 were not re-opened in this session.
- I found no fusion paper that optimises WMR weights against the exact LP worst case under pairwise constraints.

---

## Q5. Which weighted voting rules exist for K = 3 and K = 5 (enumeration)

### Takeaway
Up to relabelling, there are 2 weighted majority (decisive) rules for K = 3 and 7 for K = 5. All self-dual monotone Boolean functions with K ≤ 5 are weighted majority rules (81 labelled rules at K = 5); the first non-weighted ones appear at K = 6. So for K ≤ 5, "best weights" is a finite search over a few rules (plus signed and randomised variants), and the real content of a robust-weights paper is the structure of the answer, not the computation.

### Cited Findings
- OEIS **A001532**: "Number of NP-equivalence classes of self-dual threshold functions of n or fewer variables; number of majority (i.e., decisive and weighted) games with n players." Terms 1, 1, 2, 3, 7, 21, 135, 2470, 175428, 52980624, so a(3) = 2, a(5) = 7, a(6) = 21. References include von Neumann & Morgenstern (1944) for n ≤ 5, Gurk & Isbell (1959) for n = 6, Isbell (1959) for n = 7, and Muroga (1971) for tables up to n = 9 — [OEIS A001532](https://oeis.org/A001532)
- OEIS **A001206**: number of self-dual monotone Boolean functions of n variables: 0, 1, 2, 4, 12, 81, 2646, … (offset 0) — [OEIS A001206](https://oeis.org/A001206)
- **[OWN COMPUTATION]** Brute force over nonnegative integer weights with odd total (no ties), weights ≤ 9:
  - **K = 3:** 4 labelled WMRs (3 dictators + majority); 2 classes: (1,0,0) and (1,1,1).
  - **K = 4:** 12 labelled; 3 classes: (1,0,0,0), (1,1,1,0), (2,1,1,1).
  - **K = 5:** 81 labelled; 7 classes with minimal integer weights (1,0,0,0,0), (1,1,1,0,0), (1,1,1,1,1), (2,1,1,1,0), (2,2,1,1,1), (3,1,1,1,1) and (3,2,2,1,1).
  - **K = 6:** 1684 labelled (21 classes), against 2646 self-dual monotone functions.

  The labelled counts 4, 12, 81 equal A001206(3), A001206(4), A001206(5), so every self-dual monotone rule with K ≤ 5 is a WMR. The class counts 2, 3, 7, 21 match A001532.
- Neutral (self-dual) deterministic rules on {0,1}³ that are *not necessarily monotone*: 2⁴ = 16, of which 4 are monotone. The non-monotone ones include signed WMRs such as sign(s1 + s2 − s3) and the parity rule [OWN COMPUTATION]. Similarly, there are 2^16 neutral deterministic rules for K = 5.

### Inferences
- Under independence, the log-odds rule selects one of the 7 K = 5 classes according to which inequalities among the w_i hold. Under robust criteria, the analogous question is which of the 7 (or which signed or randomised rule) is maxmin-optimal on each region of (p, pairwise) space. That is a finite, enumerable, phase-diagram-type answer.
- Because non-monotone (negative-weight) rules can be robust-optimal (Q3 own computation), an enumeration for robust problems should include signed WMRs and randomisation, not just the 81 monotone ones.

### Gaps
- I did not verify a published statement that "all decisive simple games with ≤ 5 players are weighted". The claim rests on my computation matching OEIS counts, and is consistent with the von Neumann–Morgenstern enumeration cited in A001532.

---

## Open gaps relevant to robust weight choice

### Takeaway
"Marginals only → best single classifier" is already settled: de Oliveira, Ishii & Lin, Thm 1; Arieli et al. in the symmetric case; and the project's own LP check. It cannot be the contribution. The contribution has to come from partial dependence information. The best-supported open problem is maxmin-optimal (possibly signed or randomised) weighted voting given marginals plus pairwise joint accuracies (equivalently pairwise agreement rates), for K = 5 and general K. For K = 3 that problem is close to a special case of Balsubramani & Freund (2015) and is elementary, so it is not novel enough on its own.

### Cited Findings (evidence for each gap)
1. **Settled; do not pursue as the main result:** with marginals only, the robust-optimal "weights" put all weight on the best classifier, and no WMR improves on max p_i — [de Oliveira–Ishii–Lin](https://papers.xiao-lin.com/Robust_Information_Aggregation.pdf); [Arieli et al.](https://ar5iv.labs.arxiv.org/html/2302.03667); [Balsubramani–Freund](https://arxiv.org/pdf/1503.01811); own LP check (Q3).
2. **Open (no source found):** maxmin WMR weights under known marginals plus known pairwise joint success probabilities, for K = 5 or general K.
   - Kaniovski 2010 computes only simple-majority bounds ([preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet.pdf)).
   - Balsubramani & Freund assume the full prediction-profile distribution is known ([arXiv](https://arxiv.org/pdf/1503.01811)). That equals pairwise information only at K = 3 (own inference, Q3).
   - Questions: Is the maxmin rule always a (signed) threshold rule? When is it simple majority? How do the robust weights compare with log-odds weights? Is randomisation ever needed?
3. **Open:** sign structure. With exact marginals, can robust-optimal rules put negative weight on better-than-chance classifiers?
   - Own K = 3 computation: yes, in 22/184 random instances, with a gain of up to 0.11.
   - With lower-bound accuracy constraints, Balsubramani & Freund's σ ≥ 0 forces nonnegative weights ([arXiv](https://arxiv.org/pdf/1503.01811)).
   - A clean theorem separating "accuracy known exactly" from "accuracy lower-bounded" seems unaddressed.
4. **Open:** ambiguity sets defined by *bounds* on pairwise correlation (|ρ_ij − ρ̂_ij| ≤ ε, or ρ_ij ≤ ρ̄) for voting rules. This interpolates between the settled marginals-only case and the exact-pairwise case. Levy & Razin (2022) do the analogous thing for forecasts with a PMI bound ([IDEAS](https://ideas.repec.org/a/eee/jetheo/v199y2022ics0022053118303144.html)); Pan–Chen–Kong only treat two experts ([arXiv](https://arxiv.org/html/2311.14094v1)).
5. **Open:** dependence-invariant rules beyond homogeneous juries.
   - Kaniovski & Zaigraev's correlation-robust quota α = p(n−1)+1 exists only for exchangeable juries in a truncated Bahadur family ([preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_robust.pdf)).
   - For heterogeneous K = 3 given pairwise information, the rules with Δ(h) = 0 are invariant to the third-order moment (own inference, Q3).
   - Characterising "order-m invariant" WMRs for K = 5, and their accuracy, looks unstudied.
6. **Open:** robust weights under *structured* dependence families (opinion leader, latent groups) with *unknown* parameters. Existing work assumes known structure or studies simple or approval voting asymptotically ([Karge et al. 2024](https://iccl.inf.tu-dresden.de/w/images/e/e3/AAMAS24_OL.pdf); [Boland–Proschan–Tong 1989](https://www.cambridge.org/core/journals/journal-of-applied-probability/article/abs/modelling-dependence-in-simple-and-indirect-majority-systems/070D6335BDDDDC7AF4D70BC9B21B0B7B)). Jaffe et al. show the Bayes rule is non-linear even with a known structure ([arXiv](https://arxiv.org/pdf/1510.05830)).
7. **Partly open:** the gap between bound-minimising weights (MinCq/C-bound, tandem bound) and exact worst-case-optimal weights for small K. Bounds are known to be tight only under conditions ([Germain et al. 2015](https://ar5iv.labs.arxiv.org/html/1503.08329); [Masegosa et al. 2020](https://arxiv.org/pdf/2007.13532)).
8. **Possibly addressed but unverified:** Berg (1994) evaluates "some weighted majority decision rules under dependent voting" ([IDEAS](https://ideas.repec.org/a/eee/matsoc/v28y1994i2p71-83.html)). Shapley & Grofman (1984) "examine … interdependencies" ([Springer](https://link.springer.com/article/10.1007/BF00118940)). Nitzan & Paroush (1984) identify the optimal rule under "total interdependence" ([Springer](https://link.springer.com/article/10.1007/BF00140055)). All three must be read in full before claiming novelty on "WMRs under (known) dependence".

### Inferences
- A NeurIPS-level framing with evidence of openness: "Minimax-optimal weighted voting from pairwise agreement statistics."
  - (i) Marginals only is trivial; cite de Oliveira–Ishii–Lin.
  - (ii) With pairwise statistics, which are estimable from *unlabelled* data via agreement rates, characterise the maxmin rule.
    - For K = 3 in closed form: the feasible set is a segment in t = P(111), so the answer is a finite comparison of A(h) + t_{lo/hi}Δ(h). The connection to Balsubramani–Freund (equality constraints, so signed weights) must be acknowledged.
    - For K = 5 by LP structure: 32 cells, a 16-dimensional feasible polytope. Prove when the solution is a threshold rule, and derive the weights.
  - (iii) Compare with log-odds weights, MinCq and tandem-bound weights, and with the best single classifier.
  - (iv) Add finite-sample estimation of pairwise agreements.
- Novelty risk: K = 3 alone is elementary and close to Balsubramani & Freund. The K ≥ 5 pairwise-only game and the sign and threshold structure results are where the literature appears silent.

### Gaps
- Not checked this session: Berg (1994) full text; Shapley & Grofman's dependence section; Nitzan & Paroush (1984) definitions; Pivato (2017) weighted rules; Ladha (1995) for any optimal-weight result under correlation.
- Not searched exhaustively: the minimax-classifier / Fréchet-bound literature in operations research (Prékopa-style LP with binomial moments). Prékopa (1988) *Operations Research* 36(1):145–162 ([IDEAS](https://ideas.repec.org/a/inm/oropre/v36y1988i1p145-162.html)) and Prékopa (1990) *Operations Research* 38(2):227–239 ([INFORMS](https://pubsonline.informs.org/doi/pdf/10.1287/opre.38.2.227)) give the bounding machinery but not weight optimisation.
