# ML literature on weighted-vote guarantees under unknown dependence: second-order bounds, unlabeled agreement, weak supervision, LLM voting (novelty map for a NeurIPS paper)

Scope note: this file covers only the ML literature (NeurIPS/ICML/ICLR/COLT/AISTATS/JMLR/UAI, plus arXiv 2024–2026). The OR and voting-theory background (Matan 1996; Kuncheva et al. 2003; Kaniovski 2010; Zaigraev & Kaniovski 2012; Prékopa 1990; Boros & Prékopa 1989) is treated as already established and is not re-derived. Today's date is Oct 2026. "[unverified]" marks claims I did not confirm from the source itself. "[2024–26]" flags recent work.

---

## Q1. Second-order majority-vote bounds (C-bound, tandem, Chebyshev–Cantelli, disagreement/polarization): are they sharp, do they optimise weights, do they handle weighted votes?

### Takeaway
Every ML majority-vote bound in this line is a function of the two numbers (average error, average pairwise tandem loss or disagreement), i.e. the first two binomial moments of the error count. They all cover ρ-weighted votes, and the PAC-Bayes versions are minimised over ρ. Only the C-bound has a sharpness result, and that result holds only in a continuous-margin relaxation and only when the C-bound is at most twice the Gibbs risk. Its own authors' follow-ups state that the first-order bound beats the second-order bounds when errors are strongly correlated. That would be impossible for a bound that is sharp given both moments. I found no ML paper that uses the Boros–Prékopa/Prékopa sharp discrete-moment bounds, so idea (B) appears open. It is narrowly bracketed by a 2026 single-author arXiv paper that does sharp two-moment bounds for the exchangeable (repeated-sampling) case.

### Cited Findings
**C-bound (Lacasse et al. NIPS 2007 [venue from memory, unverified]; Germain, Lacasse, Laviolette, Marchand, Roy, JMLR 16, 2015)**
- The statement is L(MV) ≤ 1 − μ₁²/μ₂, where μ₁ and μ₂ are the first and second moments of the margin. It holds when the Gibbs risk is < 1/2. An equivalent form is 1 − (1−2R(G_Q))²/(1−2d_Q), with d_Q the expected disagreement — [Germain et al. 2015 (ar5iv)](https://ar5iv.arxiv.org/html/1503.08329); venue: [JMLR v16](https://jmlr.org/papers/v16/germain15a.html)
- **Sharpness result (Prop. 13):** "given only the first two margin moments, there exists a distribution matching them for which the C-bound equals the Bayes risk," iff 0 < μ₂ ≤ μ₁. Prop. 14 shows μ₂ ≤ μ₁ is equivalent to Gibbs risk ≤ expected disagreement, and to C-bound ≤ 2×Gibbs risk — [Germain et al. 2015](https://ar5iv.arxiv.org/html/1503.08329)
- MinCq minimises the C-bound over the posterior Q through a quadratic program, so the weights are optimised — [Germain et al. 2015](https://ar5iv.arxiv.org/html/1503.08329)
- The C-bound is hard to estimate and optimise because the second margin moment sits in the denominator. Germain et al. could minimise it only for self-complemented voter sets with aligned prior/posterior. On random forests over UCI data, the first-order bound is usually tighter than C-bounds (Lorenzen et al. 2019) — as reported by [Masegosa et al. 2020](https://ar5iv.labs.arxiv.org/html/2007.13532)

**Tandem loss bound (Masegosa, Lorenzen, Igel, Seldin, NeurIPS 2020)**
- The oracle bound is L(MV) ≤ 4·E_{ρ²}[L(h,h′)], via second-order Markov. The binary form is L(MV) ≤ 4E_ρ[L] − 2E_{ρ²}[D]. The first-order bound L(MV) ≤ 2E_ρ[L] is "tight in the worst case." When all hypotheses are identical, the second-order bound is up to 2× looser than the first-order one. The bound is minimised over ρ to reweight random-forest trees, and minimising the first-order bound instead concentrates weight and hurts test error — [Masegosa et al. (ar5iv)](https://ar5iv.labs.arxiv.org/html/2007.13532); venue: [NeurIPS 2020](https://papers.nips.cc/paper/2020/hash/386854131f58a556343e056f03626e00-Abstract.html)
- The disagreement version can use unlabeled data, but per the authors it helps only when unlabeled data is "considerably larger" than labeled data — [Masegosa et al.](https://ar5iv.labs.arxiv.org/html/2007.13532)
- The paper describes the bound itself only as "reasonably tight." It gives no sharpness theorem — [Masegosa et al.](https://ar5iv.labs.arxiv.org/html/2007.13532)

**Chebyshev–Cantelli / PAC-Bayes-Bennett (Wu, Masegosa, Lorenzen, Igel, Seldin, NeurIPS 2021)**
- The parametric oracle bound is L(MV_ρ) ≤ (E_{ρ²}[L(h,h′)] − 2μE_ρ[L] + μ²)/(½−μ)² for μ < ½. Setting μ = 0 recovers the tandem bound, and the optimal μ* recovers the C-bound. Weights are optimised with a PAC-Bayes-λ relaxation — [Wu et al. 2021 PDF](https://papers.nips.cc/paper_files/paper/2021/file/69386f6bb1dfed68692a24c8686939b9-Paper.pdf)
- **Explicit non-sharpness signal:** "in the high-error-correlation regime [E[L(h,h′)] > 0.5E[L]] the first order oracle bound is tighter than the second order oracle bounds." In addition, the empirical CCPBB bound is usually looser than TND because of union bounds over grids — [Wu et al. 2021](https://papers.nips.cc/paper_files/paper/2021/file/69386f6bb1dfed68692a24c8686939b9-Paper.pdf)
- Experiments use UCI/LibSVM datasets (SVMGuide1, Phishing, Adult, MNIST, Fashion, etc.) on 100-tree weighted random forests and heterogeneous ensembles (LDA, k-NN, trees, LR, GNB). Baselines are uniform weighting, FO, TND, CCTND and CCPBB — [Wu et al. 2021](https://papers.nips.cc/paper_files/paper/2021/file/69386f6bb1dfed68692a24c8686939b9-Paper.pdf)
- Related PAC-Bayes majority-vote work:
  - Split-kl for ternary variables: Wu & Seldin, NeurIPS 2022 — [mlanthology](https://mlanthology.org/neurips/2022/wu2022neurips-splitkl)
  - Stochastic majority votes with Dirichlet weights, giving a closed-form, differentiable risk: Zantedeschi et al., NeurIPS 2021 — [NeurIPS](https://proceedings.neurips.cc/paper/2021/hash/0415740eaa4d9decbc8da001d3fd805f-Abstract.html)
  - Diversity and generalisation for weighted MV and model averaging: Ortega, Cabañas, Masegosa, AISTATS 2022 — [PMLR](https://proceedings.mlr.press/v151/ortega22a.html)
- [2024–26] Hauptvogel & Igel (arXiv 2406.05469, 2024/25) weight deep-ensemble members by minimising the tandem PAC-Bayes bound. They report that this beats Bayes ensembles and is robust to correlated or weak members — [arXiv](https://arxiv.org/abs/2406.05469)
- [2024–26] Bastian, Leblanc, Germain, Habrard, Metzler, Morvant, Viallard (arXiv 2609.16803, Sept 2026) give disintegrated PAC-Bayes certificates for a single deterministic majority vote and a self-bounding algorithm. Tightness is characterised only by an automated review [unverified] — [Pith summary](https://pith.science/paper/2609.16803)

**Disagreement/competence bounds (Theisen, Kim, Yang, Hodgkinson, Mahoney, "When are ensembles really effective?", NeurIPS 2023 — not 2024)**
- Competence assumption: P(W_ρ ∈ [t,½)) ≥ P(W_ρ ∈ [½,1−t]) for all t ≤ ½. Under competence, EIR ≥ 0 (L(MV) ≤ E[L]), and Thm 2 gives DER ≥ EIR ≥ (2(K−1)/K)·DER − (3K−4)/K, where K is the number of classes. In the binary case this implies L(MV) ≤ 2E[L] − D (derived by the fetch tool from Thm 2; the formal corollary is in App. A.3). That is a factor-2 improvement over Masegosa's 4E[L] − 2D — [NeurIPS 2023 PDF](https://papers.neurips.cc/paper_files/paper/2023/file/30b6fa308e62ed52180c31ae3ba6bb0a-Paper-Conference.pdf); venue: [neurips.cc 2023](https://neurips.cc/virtual/2023/poster/70692)
- The paper makes no tightness claim for Thms 1–2. Weighted ρ is covered. Without competence, L(MV) = 2E[L] is attainable, so ensembling can hurt (App. C) — [Theisen et al.](https://papers.neurips.cc/paper_files/paper/2023/file/30b6fa308e62ed52180c31ae3ba6bb0a-Paper-Conference.pdf)
- Empirically, the fitted EIR–DER slopes differ from the theoretical constants (EIR ≈ DER − 1 for binary), which suggests the bounds are loose — [Theisen et al.](https://papers.neurips.cc/paper_files/paper/2023/file/30b6fa308e62ed52180c31ae3ba6bb0a-Paper-Conference.pdf)

**Polarization (Kim, Hodgkinson, Theisen, Mahoney, "How many classifiers do we need?", NeurIPS 2024)**
- η-polarization is defined by η·E[W_ρ²] ≥ P(W_ρ > ½), with η ∈ [0,4]. The bound chain is L(MV) ≤ P(W_ρ > ½) ≤ η_ρ·E_{ρ²}[L(h,h′)].
- Thm 2: L(MV) ≤ (2η(K−1)/K)(E[L] − ½E[D]).
- Corollary 1 for N weighted classifiers: L(MV) ≤ (2ηN/(N+1))(E[L] − ½E[D]).
- "Neural polarization law" (a conjecture): interpolating NN ensembles are < 4/3-polarized.
- Thm 5 on asymptotics in N: E[D_N] = (1−1/N)(D_∞ + noise). It is used to extrapolate ensemble performance from small to large N.
- Weighted votes are covered. No tightness example is given for Thms 2–4. A tightness example is shown only for the [TKY+24] bound and for 2E[L] — [arXiv html](https://arxiv.org/html/2411.00328v1); venue: [papers.cool/OpenReview (NeurIPS 2024)](https://papers.cool/venue/m5dyKArVn8@OpenReview)

### Inferences
- **The sharp bound in these bounds' own information set is unknown.** The C-bound is sharp only over margin distributions on [−1,1] (continuous). For K equal-weight voters the margin lives on a lattice, and for weighted voters it lives on the finite set of subset sums of ρ. The discrete moment problem (Prékopa 1990; Boros–Prékopa 1989) should therefore give a strictly tighter, attainable bound from exactly the same two inputs (E[#errors], E[#pairs both wrong]). Wu et al.'s "first-order beats second-order under high correlation" observation demonstrates that the published second-order bounds are not sharp given both moments.
- Kim et al.'s N/(N+1) factor shows that this line is already starting to exploit finiteness of N, but only crudely, by bounding the threshold rather than solving the moment problem.
- A sharp-bound paper could therefore (i) replace C/TND/CCPBB oracle bounds with the LP/closed-form Prékopa bound, (ii) plug it into the same PAC-Bayes concentration machinery (estimating E[L] and E[tandem] is unchanged), and (iii) minimise it over ρ. The weighted case is not a binomial-moment problem: with unequal ρ the event {MV wrong} is no longer a function of the error count. It becomes an LP over 2^K patterns, or over the support of W_ρ if only E[W_ρ] and E[W_ρ²] are used.

### Gaps
- I could not access the Lacasse et al. 2007 NIPS paper itself; the venue and year are from memory.
- I did not verify whether any C-bound paper (e.g. "On the generalization of the C-bound", arXiv 1408.1336) discusses lattice or discrete-margin sharpening. It appeared in search results but I did not read it.
- No source found states the sharp bound given (E[L], E[tandem]) for finite K. The absence of evidence relies on searches that paired "Prékopa", "Boros", "binomial moment" and "discrete moment problem" with majority vote and classifiers; none returned an ML paper.

---

## Q2. Unsupervised ensembles / accuracy estimation from agreement: identifiability, and minimax or partial-identification guarantees for weighted votes without conditional independence

### Takeaway
Two facts are classical in ML: the pairwise agreement identity, and the underdetermination of higher-order dependence (Platanios et al., UAI 2014). The spectral/tensor line (Parisi 2014; Jaffe 2015/16; Zhang 2016) obtains point identification only by assuming conditional independence, or a known/sparse dependence structure, plus a sign/"better than random" assumption to break the ±1 label-flip ambiguity. The closest worst-case results without CI are two. Balsubramani & Freund (COLT 2015; NeurIPS 2016) give a transductive minimax aggregator that uses per-classifier accuracy lower bounds plus the full unlabeled prediction matrix; the optimum is a clipped weighted vote. Mazzetto et al. (AISTATS 2021) give an LP worst-case majority-vote error from labeled error rates plus unlabeled pairwise disagreements, but for unweighted votes over selected subsets. Nobody states the even/odd (label-flip-symmetric) decomposition explicitly, or that self-dual weighted-vote accuracy depends only on the odd part. Nobody gives sharp identified sets or minimax deterministic weights for weighted votes from (pᵢ, agreement).

### Cited Findings
**Agreement identities and underdetermination**
- Platanios, Blum, Mitchell (UAI 2014) give the pairwise identity a_{ij} = 1 − e_i − e_j + 2e_{ij}. The subset version, via inclusion–exclusion, is a_A = e_A + 1 + Σ_k (−1)^k Σ_{|I|=k} e_I. N functions give 2^N − N − 1 agreement equations against 2^N − 1 unknown error rates, so the system is "underdetermined." Exact identification requires independent errors with N ≥ 3 and e_i < 0.5. Otherwise they minimise a dependence penalty Σ_A (e_A − Π e_i)² subject to the agreement equations — [UAI 2014 PDF](https://www.auai.org/uai2014/proceedings/individuals/313.pdf)
  - Note the overlap with idea (A): r_ij = P(Z_i = Z_j = 1) = 1 − e_i − e_j + e_ij, so (p_i, p_j, a_ij) identify r_ij. This is a direct re-expression of Platanios' identity. The novelty of (A) cannot rest on this identity.
- Follow-ups: Platanios, Dubey, Mitchell, "Estimating accuracy from unlabeled data: a Bayesian approach," ICML 2016 — [mlanthology](https://mlanthology.org/icml/2016/platanios2016icml-estimating); Platanios, Poon, Mitchell, Horvitz, "…a probabilistic logic approach," NeurIPS 2017 — [NeurIPS](https://proceedings.neurips.cc/paper/2017/hash/95f8d9901ca8878e291552f001f67692-Abstract.html) [details unverified]

**Spectral / tensor unsupervised ensembles (CI and sign assumptions)**
- Parisi, Strino, Nadler, Kluger, "Ranking and combining multiple predictors without labeled data" (PNAS 2014; venue unverified from page). They assume conditional independence, under which the off-diagonal covariance is rank one. The leading eigenvector ∝ balanced accuracies and is identified only "up to a sign ambiguity." The spectral meta-learner weights classifiers by the eigenvector entries. A "cartel" of colluding classifiers gives rank two; SML is robust to a small cartel — [arXiv/ar5iv](https://ar5iv.labs.arxiv.org/html/1303.3257)
- Jaffe, Nadler, Kluger, AISTATS 2015 assume CI and that most classifiers have balanced accuracy > 0.5. In their words, "Assumption (iii) or a variant thereof is needed, given an inherent ±1 sign ambiguity." The off-diagonal covariance matrix and third-order tensor are rank one, which identifies class imbalance plus sensitivities and specificities — [ar5iv](https://ar5iv.arxiv.org/html/1407.7644); [PMLR v38](https://proceedings.mlr.press/v38/jaffe15.html)
- Jaffe, Fetaya, Nadler, Jiang, Kluger, AISTATS 2016 ("Unsupervised ensemble learning with dependent classifiers"). They relax CI with a model of dependent classifiers, detect strongly dependent ones, and build an improved meta-learner. The abstract mentions no worst-case guarantees. I could not verify the exact model (I believe it is groups or a latent tree) — [PMLR v51](https://proceedings.mlr.press/v51/jaffe16.html)
- Zhang, Chen, Zhou, Jordan, "Spectral methods meet EM," JMLR 17 (2016): spectral initialisation followed by EM for Dawid–Skene. It "achieves the optimal convergence rate up to a logarithmic factor" (the CI assumption is from Dawid–Skene; I did not confirm it from that page) — [JMLR](https://www.jmlr.org/papers/v17/14-511.html)
- Berend & Kontorovich, "Consistency of weighted majority votes," NeurIPS 2014. They assume independent experts. With known competences the optimal rule is the Nitzan–Paroush log-odds weighting, with error ≤ exp(−Φ/2) and a matching-order lower bound. With unknown competences, no frequentist-optimal empirical rule exists — [ar5iv](https://ar5iv.arxiv.org/html/1312.0451); [NeurIPS 2014](https://papers.nips.cc/paper_files/paper/2014/hash/3cf61323ed04cc4581f435dab1256ab8-Abstract.html)
- [2024–26] Corrada-Emmanuel, "Algebraic Evaluation Theorems" (arXiv 2412.16238, 2024/25). Under error independence, agreement-pattern counts of ≥3 binary classifiers give exactly two solutions for the correctness statistics. This is a label-flip-type two-fold ambiguity, resolved by picking the solution with the higher total accuracy. Correlated errors generally make the closed form irrational, which serves as an "alarm" — [arXiv html](https://arxiv.org/html/2412.16238v2)

**Worst-case / minimax aggregation without CI (closest prior art to idea A)**
- Balsubramani & Freund, "Optimally combining classifiers using unlabeled data," COLT 2015 (PMLR v40). This is transductive. The adversary chooses test labels z ∈ [−1,1]^n subject only to (1/n)Fz ≥ b, where F is the full unlabeled prediction matrix and b holds per-classifier correlation (accuracy) lower bounds from labeled data. There is **no independence assumption**. The minimax prediction is a **clipped weighted vote** g = clip(xᵀσ*), with σ* ≥ 0 minimising a convex slack function γ(σ) = (1/n)Σ[|xⱼᵀσ| − 1]₊ − bᵀσ. There are no explicit pairwise constraints; dependence enters through F — [arXiv PDF v2](https://arxiv.org/pdf/1503.01811v2); [PMLR](https://proceedings.mlr.press/v40/Balsubramani15.html)
- Balsubramani & Freund, "Optimal binary classifier aggregation for general losses," NeurIPS 2016. The minimax optimal prediction for general monotone losses is a sigmoid-like link of a weighted vote, again with no independence assumption — [ar5iv](https://ar5iv.labs.arxiv.org/html/1510.00452); [NeurIPS 2016](https://proceedings.neurips.cc/paper/2016/hash/eaa52f3366768bca401dca9ea5b181dd-Abstract.html)
- Balsubramani & Freund, "Scalable semi-supervised aggregation of classifiers" (NeurIPS 2015 per corpus path; unverified). HedgeClipper handles specialists and abstention, is optimised by SGD, and is compared against RF, AdaBoost and LR on UCI/LibSVM/Kaggle (e.g. SUSY) — [text](https://object.cloud.sdsc.edu/v1/AUTH_da4962d3368042ac8337e2dfdd3e7bf3/ml-papers-txt/NeurIPS/2015/scalable_semisupervised_aggregation_of_classifiers__60ec1410.txt)
- Mazzetto, Sam, Park, Upfal, Bach, "Semi-supervised aggregation of dependent weak supervision sources with performance guarantees," AISTATS 2021 (PMLR v130):
  - Inputs are individual error rates ε_i from labeled data and pairwise disagreements d_ij from unlabeled data (or the full output distribution Q).
  - **No assumption on the joint distribution.** The worst-case MV error of a subset is computed by an LP over {0,1}^k correctness patterns.
  - Prop. 5 gives a closed form for 3 labelers: ε_i + ε_j + ε_k − ΣD/2 − 2q. In the equal case it is ε + d/2 if d ≤ ε, else 3(ε − d/2). The worst case is ≤ ε iff d ≥ (4/3)ε.
  - Prop. 2: with dependent labelers, MV over all n can have error ≥ the median individual error.
  - Prop. 3: **given only error rates, no function beats the best single labeler in the worst case.**
  - Unweighted MV over greedily selected subsets only; weighted voting appears only as background (Nitzan–Paroush).
  - Experiments: 45 binary tasks from AwA2. Baselines are MV, MV with flips, semi-supervised Dawid–Skene and ALL.
  - It does not cite Kuncheva or Matan — [NSF PAR full text](https://par.nsf.gov//servlets/purl/10273522); [mlanthology](https://mlanthology.org/aistats/2021/mazzetto2021aistats-semisupervised)

**Disagreement as a proxy for error (not identification)**
- Jiang, Nagarajan, Baek, Kolter, "Assessing generalization of SGD via disagreement" (ICLR 2022): test error is estimated by disagreement between two SGD runs on unlabeled data, explained via ensemble calibration — [arXiv](https://arxiv.org/abs/2106.13799); venue: [mlanthology ICLR 2022](https://mlanthology.org/iclr/2022/jiang2022iclr-assessing)
- Baek, Jiang, Raghunathan, Kolter, "Agreement-on-the-line" (NeurIPS 2022): ID-vs-OOD agreement is linear with the same slope and bias as accuracy, which allows OOD accuracy estimation from unlabeled data. It is observed only for NNs — [arXiv](https://arxiv.org/abs/2206.13089); venue: [mlanthology NeurIPS 2022](https://mlanthology.org/neurips/2022/baek2022neurips-agreementontheline)

**Partial identification framing (closest framing precedent)**
- [2024–26] Maia Polo, Maity, Yurochkin, Banerjee, Sun, "Weak supervision performance evaluation via partial identification," NeurIPS 2024. Fréchet bounds on accuracy/F1 of a classifier, given P(X,Z) (unlabeled) and P(Y,Z) (label model). The bounds are sharp (attained), computed through an OT dual, with asymptotic CIs. It does **not** optimise aggregation weights. Experiments are on WRENCH — [arXiv html](https://arxiv.org/html/2312.04601v2); venue: [NeurIPS 2024](https://proceedings.neurips.cc/paper_files/paper/2024/hash/f4c6bec746b0aeca8c2cd15096f1ad1f-Abstract.html)
- [2024–26] Chen, Rambachan, Tamer, "Partial identification from LLM prompts" (econ.EM, arXiv 2606.15031, June 2026). Prevalence θ with LLM errors **arbitrarily dependent** given the true label. Without identifying restrictions the identified set is [0,1]; identifying content must come from externally calibrated accuracy-type inputs. They list "weighted ensembles" among nested calibrated scores but derive no bounds for weighted voting — [arXiv html](https://arxiv.org/html/2606.15031)

### Inferences
- **Idea (A), the observational-equivalence part, is half-known.** The sign ambiguity is acknowledged throughout the spectral literature, and Platanios shows underdetermination. But the specific structural statement is not stated anywhere I found: with ±1 coding, Z_iZ_j = f_if_j, so all even-order correctness moments are identified from unlabeled predictions, while odd-order ones carry the unobservable factor Y. Combined with "the accuracy of a self-dual weighted vote is a functional of the odd part only," this is a clean, citable theorem-level contribution.
- **Idea (A), the "minimax aggregation" part, is substantially anticipated** by Balsubramani & Freund 2015. Their feasible set (observed f-law plus first-order accuracy constraints) is exactly the identified set when the full unlabeled prediction law is known. Their aggregator is randomized/soft (clipped), and it is transductive.
- What is new in (A): deterministic self-dual weighted votes; population (inductive) identified sets; the use of only pairwise agreement (which is what can be estimated reliably for large K); sharp intervals [L(w), U(w)] for any given weight vector; minimax weights argmax_w L(w); the gap between the minimax deterministic vote and the B&F randomized optimum; and finite-sample confidence sets.
- Mazzetto Prop. 3 (and de Oliveira–Ishii–Lin, see Q5) already establishes that with marginals only the robust rule is a dictator. The interesting regime for (A) is therefore marginals plus pairwise information, which Mazzetto treats only for unweighted subset MVs. Reviewers who know Mazzetto 2021 will ask exactly how far (A) goes beyond it. The answer must be weighted votes, sharpness/closed forms, the even/odd theorem, and minimax weights.

### Gaps
- Not verified: the exact dependence model in Jaffe et al. 2016; the details of Platanios 2016/2017; and whether any 2025–26 paper (beyond those found) computes LP-based worst-case accuracy for weighted votes. Several targeted searches found none.
- I did not check Corrada-Emmanuel's NTQR papers (arXiv 2312.05392, 2509.08593, 2510.00821) for an explicit even/odd statement. These "logical consistency/no-knowledge alarm" works are the most likely place for overlapping algebra and should be read before submission.

---

## Q3. Weak supervision with dependent labeling functions: assumptions and worst-case guarantees

### Takeaway
The Snorkel line models dependence parametrically and learns it from unlabeled data:
- Ratner et al. NeurIPS 2016: data programming.
- Bach et al. ICML 2017: structure learning via ℓ₁ pseudolikelihood.
- Varma et al. ICML 2019: robust-PCA structure.
- Ratner et al. AAAI 2019: MeTaL.
- Fu et al. ICML 2020: triplet/FlyingSquid on a user-provided or learned Ising graph.

Their guarantees are consistency and sample-complexity results under a correctly specified (or KL-bounded misspecified) model, not worst-case bounds. Worst-case guarantees come from a separate "adversarial label learning" line (Arachie & Huang; Mazzetto et al. 2021 ×2), which uses per-source constraints, and at most pairwise disagreements, with unweighted aggregation or learned classifiers.

### Cited Findings
- Ratner, De Sa, Wu, Selsam, Ré, "Data programming," NeurIPS 2016. A generative model of labeling functions, with parameter recovery "in a few settings" — [arXiv](https://arxiv.org/abs/1605.07723v3); venue: [neurips.cc 2016](https://neurips.cc/virtual/2016/poster/7031)
- Bach, He, Ratner, Ré, "Learning the structure of generative models without labeled data," ICML 2017. ℓ₁-regularised marginal pseudolikelihood. The unlabeled data needed "scales sublinearly in the number of possible dependencies" — [PMLR v70](https://proceedings.mlr.press/v70/bach17a.html)
- Varma, Sala, He, Ratner, Ré, "Learning dependency structures for weak supervision models," ICML 2019. Robust-PCA-based, assuming sparsity of the dependency structure. Sample complexity is sublinear or logarithmic in m, with an information-theoretic lower bound and no labels — [PMLR v97](https://proceedings.mlr.press/v97/varma19a.html)
- Ratner et al., "Training complex models with multi-task weak supervision," AAAI 2019 (MeTaL; matrix completion given the dependency graph) [details unverified] — [AAAI](https://ojs.aaai.org/index.php/AAAI/article/view/4403)
- Fu, Chen, Sala, Hooper, Fatahalian, Ré, "Fast and Three-rious" (FlyingSquid), ICML 2020:
  - Binary Ising model over (Y, λ), with a user-provided or estimated dependency graph that is triangulated.
  - Closed-form triplet estimation.
  - The analysis does not assume the latent model is exactly correct, and a KL misspecification analysis is included. No worst-case guarantee for a misspecified graph was found in the excerpt — [ar5iv](https://ar5iv.labs.arxiv.org/html/2002.11955); venue: [ICML 2020](https://icml.cc/virtual/2020/poster/53479)
- Arachie & Huang, "Adversarial label learning" (AAAI 2019 [venue unverified]); "A general framework for adversarial label learning," JMLR 22 (2021). The adversary picks labels within expert-provided error bounds of the weak signals. The primal objective upper-bounds the true error if the bounds hold. The authors say it is "robust against dependent and biased errors"; correlated copies of a signal leave the optimum unchanged — [ar5iv](https://ar5iv.labs.arxiv.org/html/1805.08877); [JMLR](https://www.jmlr.org/beta/papers/v22/20-537.html)
- Mazzetto, Cousins, Sam, Bach, Upfal, "Adversarial multi class learning under weak supervision with performance guarantees," ICML 2021 (PMLR v139):
  - Constraints are per-source only (risks or confusion-type statistics); **no pairwise statistics**. Sources may be "arbitrarily correlated."
  - The guarantee is risk ≤ best-in-class + (feasible-set diameter) + complexity terms.
  - The union-bound analysis is acknowledged as loose for correlated sources — [NSF PAR](https://par.nsf.gov//servlets/purl/10274517); [PMLR](https://proceedings.mlr.press/v139/mazzetto21a.html)
- Mazzetto et al. AISTATS 2021 is the closest; see Q2 for details — [NSF PAR](https://par.nsf.gov//servlets/purl/10273522)
- WRENCH (NeurIPS 2021 Datasets & Benchmarks) is the standard weak-supervision benchmark. Maia Polo et al. 2024 evaluate on it — [NeurIPS D&B](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/hash/1c9ac0159c94d8d0cbedc973445af2da-Abstract-round2.html)

### Inferences
- No weak-supervision paper I found provides sharp worst-case accuracy for a weighted vote given source accuracies plus pairwise agreements. The adversarial line uses first-order constraints, or pairwise constraints only for unweighted subset MVs.
- The structure-learning line (Bach, Varma, FlyingSquid) needs sparsity or a correct graph. Idea (A) could be pitched as the assumption-free complement: what is guaranteed if the graph is wrong.
- WRENCH is a natural benchmark where many labeling functions are dependent. However, labeling functions abstain, so the theory would need an abstention extension (cf. B&F specialists).

### Gaps
- I did not verify the theory sections of Ratner 2016/2019 (exact identifiability conditions) or Varma 2019 (exact recovery conditions).

---

## Q4. LLM-era relevance: non-monotone voting, correlated LLM errors, juries/judge ensembles — is there a theory gap that robust weighted voting fills?

### Takeaway
Since 2024 there has been a flood of empirical evidence that LLM errors are strongly correlated across samples, models and judges. Theory has started to appear in 2025–26, but it splits into two kinds:
- (i) Exchangeable repeated sampling from one protocol: Chen et al. NeurIPS 2024; Liu 2026 ×2, which include **sharp two-moment bounds and non-monotone vote curves**; Bay & Yearick 2026.
- (ii) Parametric dependence models for multi-model or judge aggregation: Ising (Balasubramanian et al. 2026), confounder low-rank (CARE 2026), Gaussian copula (Turkmen et al. 2026), and conditional independence (Ai et al. ICML 2026).

I found no paper giving sharp, assumption-free guarantees or minimax weights for **weighted** votes over **heterogeneous** LLMs or judges with arbitrary dependence. That is the gap ideas (A) and (B) fill. Idea (C) is now heavily covered for the exchangeable case, but not for heterogeneous worst-case ensembles under pairwise information.

### Cited Findings
**Non-monotone / plateauing voting**
- Chen, Davis, Hanin, Bailis, Stoica, Zaharia, Zou, "Are more LLM calls all you need? Towards the scaling properties of compound AI systems," NeurIPS 2024. Vote and Filter-Vote accuracy "can first increase but then decrease" with the number of calls. More calls help easy queries (per-call accuracy > ½) and hurt hard ones. The analysis assumes calls are **i.i.d. given the query**, uses an unweighted vote and does not discuss correlated calls — [neurips.cc](https://neurips.cc/virtual/2024/poster/93777); [ar5iv](https://ar5iv.labs.arxiv.org/html/2403.02419)
- Brown et al., "Large language monkeys" (arXiv 2407.21787, 2024). Coverage scales log-linearly, but majority voting and reward models "plateau beyond several hundred samples" — [arXiv](https://arxiv.org/abs/2407.21787v3)
- [2024–26] **Liu (York U.), "Two calls, two moments, and the vote-accuracy curve of repeated LLM inference"** (arXiv 2605.03379, v3 Oct 2026). This is the closest recent work to (B) and (C):
  - **Sharp** finite-budget bounds on majority-vote accuracy given μ = E[Q] and ν = E[Q²] (two-call both-correct rate). They are computed via quadratic minorants/majorants of the binomial-majority polynomial P_n(q), and extremes are attained by ≤3-point laws.
  - Thm 4.2 gives a closed form for 3 votes. Gains from 3 votes are guaranteed iff ρ < 1 − 1/(2μ). Infinite-vote sharp bounds are also given.
  - Assumes **conditional independence of calls given a per-example Q** (de Finetti mixture). **Unweighted, single protocol; heterogeneous models and weighted votes are out of scope.**
  - Cites Karlin–Studden, Bertsimas–Popescu, Winkler and Pinelis, but **not Prékopa, Kuncheva/Matan or the PAC-Bayes tandem/C-bound line**.
  - Experiments: QNLI/QQP with Llama-3.1-8B, Phi-4-mini and Qwen2.5-7B; a GSM8K supplement — [arXiv PDF](https://arxiv.org/pdf/2605.03379)
- [2024–26] Liu, "When can voting help, hurt, or change course? Exact structure of binary test-time aggregation" (arXiv 2605.05592, 2026):
  - The vote curve determines exactly a signed "voting signature," not the latent law.
  - Example 4.1 is a law whose curve **changes direction infinitely often**, so no general eventual-monotonicity theorem exists.
  - Infinite exchangeability; weighted votes "explicitly outside scope"; no experiments — [arXiv PDF](https://arxiv.org/pdf/2605.05592)
- [2024–26] Bay & Yearick, "When more sampling hurts: the modal ceiling and correlation ceiling of test-time scaling" (arXiv 2606.28661, June 2026):
  - Design effect n_eff = n/[1+(n−1)ρ], with ceiling 1/ρ.
  - Self-consistency converges to the modal-hit rate, and per-problem accuracy → 0 when the mode is wrong.
  - Benchmark-level curves plateau rather than decline. There are no worst-case bounds and no weighted votes — [alphaxiv](https://www.alphaxiv.org/abs/2606.28661.md)
- Other ensemble-size results:
  - Aeeneh, Zlatanov, Yu (arXiv 2309.09564, 2023): multi-class MV with independent voters; error decays or **grows exponentially** in the number of voters depending on conditions — [arXiv](https://arxiv.org/abs/2309.09564)
  - Bax (arXiv 2104.11833, 2021): "any number of voters may minimize the error rate" for random-subset voting — [arXiv](https://arxiv.org/abs/2104.11833)

**Correlated LLM errors (empirical)**
- Kim, Garg, Peng, Garg, "Correlated errors in large language models," ICML 2025:
  - 350+ LLMs studied. When two models are both wrong, they agree on the same wrong answer 60% of the time on one leaderboard dataset.
  - Shared provider or architecture increases correlation, and larger, more accurate models are more correlated.
  - Effects shown on LLM-as-judge and on hiring (algorithmic monoculture) — [PMLR v267](https://proceedings.mlr.press/v267/kim25e.html)
- Goel et al., "Great models think alike and this undermines AI oversight," ICML 2025. The CAPA similarity metric; judges favour similar models; mistakes become more alike as capability grows — [PMLR v267](https://proceedings.mlr.press/v267/goel25b.html)
- [2024–26] Kohli (Apple), "Nine judges, two effective votes" (arXiv 2605.29800, May 2026):
  - A 9-judge, 7-family panel has n_eff ≈ 2.18 (mean pairwise error φ = 0.391).
  - The panel reaches 72.0% on ChaosNLI-MNLI, versus 94% predicted under independence, and the best single judge sometimes beats the panel.
  - Dawid–Skene and **oracle accuracy-weighted voting** close at most 11% of the "Condorcet gap."
  - Benchmarks: ChaosNLI (MNLI/SNLI/αNLI) and RewardBench — [arXiv html](https://arxiv.org/html/2605.29800v1)
- [2024–26] Denisov-Blanch, Kazdan, …, Koyejo, "Consensus is not verification" (arXiv 2603.06612 [date unverified]):
  - Majority, confidence-weighted and surprisingly-popular aggregation show no consistent gains over single samples on truthfulness tasks at ~25× cost.
  - Models agree above chance even on random strings.
  - "No rule based only on internal signals can reliably separate a unanimous correct answer from a unanimous wrong one" — [arXiv html](https://arxiv.org/html/2603.06612)
- [2024–26] J. Chen, "When does combining language models help? A co-failure ceiling…" (arXiv 2606.27288, June 2026):
  - Accuracy of any member-selecting policy is ≤ 1 − β, where β = P(all models wrong).
  - **Pairwise correlation cannot identify β**, since equal marginals and pairwise statistics can give different all-wrong rates.
  - A Gaussian copula underprices β about 2.5×. The study covers 67 models — [arXiv](https://arxiv.org/abs/2606.27288)
- [2024–26] Kim (Aidentyx), "Are diversity metrics measuring diversity?" (arXiv 2607.20768, July 2026). Over 31,900 LLM subsets on MMLU-Pro and TruthfulQA, simple voting beats the best member in only 9.98% of size-3 subsets. Diversity metrics are nearly collinear with mean accuracy, and the only stable signal is pairwise co-failure — [arXiv html](https://arxiv.org/html/2607.20768v1)
- [2024–26] Nitarach, AIMO 3 report (arXiv 2603.27844, Mar 2026). Uses N_eff = N/(1+(N−1)ρ) and finds that prompt-diversity interventions failed — [arXiv html](https://arxiv.org/html/2603.27844v1)

**LLM juries and aggregation methods**
- Verga et al. (Cohere), "Replacing judges with juries" (PoLL; arXiv 2404.18796, 2024; venue unverified):
  - Three diverse small judges with max-vote or average pooling beat a GPT-4 judge on κ (e.g. NQ 0.763 vs 0.627) at about 7× lower cost.
  - Voting is unweighted, and there is **no theory of judge correlation** — [arXiv html](https://arxiv.org/html/2404.18796v2)
- [2024–26] Ai, Pan, Simchi-Levi, Tambe, Xu, "Beyond majority voting: LLM aggregation by leveraging higher-order information," **ICML 2026** (PMLR 306):
  - Optimal Weight (OW) uses first-order accuracies. Inverse Surprising Popularity (ISP) uses second-order conditional answer probabilities P(A_i|A_j).
  - Assumes **conditional independence given the true label** (with a difficulty extension) and random answer shuffling.
  - **No minimax guarantee** beyond a shuffling argument.
  - Benchmarks: UltraFeedback, MMLU, ARMMAN — [arXiv PDF](https://arxiv.org/pdf/2510.01499); venue: [PMLR v306](https://proceedings.mlr.press/v306/ai26d.html)
- [2024–26] Balasubramanian, Podkopaev, Kasiviswanathan (UC Davis/AWS), "Dependence-aware label aggregation for LLM-as-a-judge via Ising models" (arXiv 2601.22336, Jan 2026):
  - With class-independent couplings, the Bayes rule is linear in votes with correlation-adjusted weights; with class-dependent couplings it is quadratic.
  - CI methods can flip the Bayes label even with matched marginals. Inference is by EM with pseudolikelihood.
  - **No worst-case bounds.** Tasks: relevance, toxicity and summarisation with six commercial judges; 0.90 vs 0.82 for CI-weighted MV on relevance — [alphaxiv](https://www.alphaxiv.org/abs/2601.22336)
- [2024–26] Zhao, Shin, Huang, Namburi, Sala, "CARE: Confounder-aware aggregation for reliable LLM evaluation" (arXiv 2603.00039, 2026). Sparse-plus-low-rank latent confounders plus tensor decomposition; identifiability up to a symmetry-breaking step (anchors); 12 benchmarks, with up to 26.8% error reduction versus MV and averaging — [arXiv html](https://arxiv.org/html/2603.00039)
- [2024–26] Turkmen, Buyukates, Bastopcu, "Don't always pick the highest-performing model" (arXiv 2602.08003, Feb 2026):
  - Under independence top-k is optimal; under Gaussian-copula correlation the MAP error has a positive floor.
  - Greedy mutual-information selection. MAP aggregation, no weighted vote and no worst-case guarantee.
  - Benchmarks: MEDMCQA, MMLU, IMDB — [arXiv PDF](https://arxiv.org/pdf/2602.08003)
- [2024–26] Acharya, Pan, Verkhovsky (AWS), "ROPOLL: Robust panel of LLM judges" (arXiv 2606.30931, June 2026). Huber ε-contamination robustness via the geometric median, assuming conditionally independent judges with an equicorrelated extension. Benchmarks: HelpSteer 2/3 and UltraFeedback — [arXiv PDF](https://arxiv.org/pdf/2606.30931)

### Inferences
- **The theory gap is specific and real.** Every 2024–26 multi-model aggregation theory either assumes a dependence model (CI, Ising, Gaussian copula, low-rank confounders, equicorrelation) or restricts to exchangeable repeated calls of one model (Liu; Chen et al.; Bay & Yearick).
- Two recent results directly motivate a robust approach:
  - Kohli shows that oracle accuracy-weighted voting fails under correlation.
  - J. Chen shows pairwise statistics cannot pin down higher-order co-failure, so a model-based fill-in is arbitrary.
- (A)+(B) offers a principled replacement: sharp intervals for accuracy, minimax weights and certified "add/drop this model" decisions, using exactly the observables practitioners have (small labeled set for accuracies, large unlabeled pool for agreement).
- For (C), worst-case non-monotonicity is no longer surprising as a phenomenon; Chen et al., Liu and the marginals-only Kuncheva result all cover versions of it. What remains new is a worst-case statement over heterogeneous members given pairwise information, e.g. "adding a model with accuracy p and agreement a_j can strictly lower the sharp lower bound." Without a weighted or heterogeneous angle, (C) alone would likely be judged incremental relative to Liu 2605.03379/2605.05592 and Chen et al. 2024.

### Gaps
- Venue status of the Liu papers, Kohli, J. Chen, Bay & Yearick, CARE, ROPOLL and Balasubramanian et al. is unknown (arXiv only as of my reading). Some may be under review at NeurIPS/ICLR 2026–27.
- I could not access the full statements of J. Chen's finite-sample certificate or Kohli's "Condorcet null" construction.

---

## Q5. Recent (2023–2026) "distributionally robust / minimax ensemble weighting / robust aggregation with unknown dependence"

### Takeaway
In ML proper there is almost nothing under the name "distributionally robust ensemble weighting" that treats classifier dependence as the ambiguity set. Search hits for "minimax majority vote" turned out to concern dependence among training samples (Shihab et al. 2026) or adversarial judge contamination (ROPOLL). The closest conceptual literature is in economics and computation on robust information/forecast aggregation (Arieli et al. 2018; Levy & Razin JET 2022; de Oliveira–Ishii–Lin; Pan–Chen–Kong WWW 2024; Guo et al.; Chen–Peng–Tang 2026). It establishes a key negative baseline: with binary states and actions and known marginals only, the robust optimum relies on a single source.

### Cited Findings
- de Oliveira, Ishii, Lin, "Robust aggregation of correlated information" (arXiv 2106.00088, v2 Sept 2024). The decision maker knows each source's marginal information structure but not the correlation and maximises the worst case. **With binary states and actions, a robustly optimal strategy relies on a single source.** For general problems, the number of sources needed is bounded — [arXiv](https://arxiv.org/abs/2106.00088)
- Levy & Razin, "Combining forecasts in the presence of ambiguity over correlation structures," JET 199 (2022). The correlation ambiguity set is bounded through pointwise mutual information, and the set of predictions is determined by the naive-Bayes forecast plus the correlation bound — [RePEc](https://ideas.repec.org/a/eee/jetheo/v199y2022ics0022053118303144.html)
- Pan, Chen, Kong, "Robust decision aggregation with second-order information" (WWW 2024 per DOI 10.1145/3589334.3645384 [venue from search listing]). Two experts, binary state and action. Over all information structures, second-order information (each expert's prediction of the other) gives no robustness benefit, with regret 0.5 achieved by "follow the first expert." It helps in some conditionally independent families — [arXiv html](https://arxiv.org/html/2311.14094v1)
- Guo, Hartline, Huang, Kong, Shah, Yu, "Algorithmic robust forecast aggregation" (arXiv 2401.17743; venue unverified). FPTAS for minimax-regret aggregators over families of information structures — [arXiv](https://arxiv.org/pdf/2401.17743)
- [2024–26] Chen, Peng, Tang, "Prior-agnostic robust forecast aggregation" (arXiv 2604.24517, May 2026). Over general (arbitrarily correlated) structures the minimax regret is 1/4, achieved by outputting ½, which is trivial. Positive results hold only for CI or Blackwell-ordered classes — [arXiv](https://arxiv.org/pdf/2604.24517)
- [2024–26] Huang, Xia, Kong, "Decision aggregation under quantal response" (arXiv 2603.13807 [date unverified]). Majority voting is minimax-regret optimal when rationality is low, under conditionally i.i.d. signals. Correlated signals are named as future work — [arXiv html](https://arxiv.org/html/2603.13807v1)
- [2024–26] Shihab, Akter, Sharma, "Minimax optimality and spectral routing for majority-vote ensembles under Markov dependence" (arXiv 2604.13414, Apr 2026). The dependence studied is among **training samples** (Markov chains), not among classifiers; not relevant except as a name collision — [arXiv](https://arxiv.org/abs/2604.13414)
- [2024–26] Bektas & Can (arXiv 2511.21465, Nov 2025). Ensemble size versus the probability of linear independence of votes, including weighted MV; geometric and heuristic, with no worst-case accuracy — [arXiv](https://arxiv.org/abs/2511.21465)

### Inferences
- The economics results are the right "impossibility" anchor:
  - With only marginal accuracies, minimax weights put all weight on one classifier (de Oliveira et al.; Mazzetto Prop. 3). This also matches the known marginals-only voting bounds.
  - A paper on (A) should show that adding identified pairwise information (agreement) breaks this dictatorship, and characterise when. To my knowledge no source does this for weighted votes with K ≥ 3.
  - Pan–Chen–Kong's "second-order information does not help over all structures" uses a different notion of second order (beliefs about others' reports, two experts). It must be distinguished explicitly to pre-empt reviewer confusion.
- The phrase "distributionally robust weighting" is essentially unused for this problem in ML venues, which is a positioning opportunity. It also means reviewers from the DRO community may expect Wasserstein or φ-divergence ambiguity sets. The paper should explain why the moment/agreement-defined ambiguity set (a partial-identification set) is the natural one.

### Gaps
- I did not search SSRN/econ journals exhaustively; other researchers may cover voting-theory and econ work (e.g. Kaniovski–Zaigraev optimal jury weights with correlation; Peleg–Zamir dependent Condorcet).
- I found no ML paper titled or framed as "minimax weights for weighted majority vote under unknown classifier dependence." This is absence of evidence from roughly 10 targeted queries, not proof.

---

## Gaps and positioning for a NeurIPS submission

### Takeaway
The strongest novel package is (A)+(B) for weighted votes over heterogeneous classifiers or LLMs:
- (i) the even/odd identification theorem;
- (ii) sharp identified intervals for the accuracy of any weighted vote from (labeled p_i, unlabeled agreements);
- (iii) minimax (max–min) weights, with a "dictator → non-dictator" phase transition;
- (iv) sharp discrete-moment (Prékopa) replacements for the C-bound/tandem/CCPBB bounds, with a demonstration that the published bounds are not sharp.

(C) should be a corollary or showcase (worst-case non-monotonicity in K under pairwise information for heterogeneous members), not the headline, because Chen et al. (NeurIPS 2024) and Liu (arXiv 2026, ×2) cover non-monotone voting for exchangeable calls, and Kuncheva 2003 covers marginals-only.

### Cited Findings
**Already done (must cite and differentiate)**
- The pairwise identity linking agreement and joint error (therefore r_ij identification) is Platanios et al. UAI 2014 — [PDF](https://www.auai.org/uai2014/proceedings/individuals/313.pdf). Underdetermination of higher-order terms is shown in the same paper.
- The sign/label-flip ambiguity of unlabeled methods is acknowledged by Jaffe et al. 2015 — [ar5iv](https://ar5iv.arxiv.org/html/1407.7644), Parisi et al. 2014 — [ar5iv](https://ar5iv.labs.arxiv.org/html/1303.3257), and Corrada-Emmanuel 2024/25 — [arXiv](https://arxiv.org/html/2412.16238v2)
- Transductive minimax aggregation with accuracy constraints, no independence assumption and a clipped weighted vote: Balsubramani & Freund COLT 2015 — [arXiv](https://arxiv.org/pdf/1503.01811v2); NeurIPS 2016 — [ar5iv](https://ar5iv.labs.arxiv.org/html/1510.00452)
- LP worst-case unweighted MV error from labeled ε_i plus unlabeled pairwise disagreement, with a 3-labeler closed form and "no function beats the best single labeler given only error rates": Mazzetto et al. AISTATS 2021 — [NSF PAR](https://par.nsf.gov//servlets/purl/10273522)
- Sharp two-moment MV bounds and possible negative vote gains for exchangeable repeated LLM calls: Liu 2026 — [arXiv](https://arxiv.org/pdf/2605.03379); infinitely-oscillating vote curves: Liu 2026 — [arXiv](https://arxiv.org/pdf/2605.05592)
- Non-monotone Vote with call count: Chen et al. NeurIPS 2024 — [neurips.cc](https://neurips.cc/virtual/2024/poster/93777)
- Robust single-source optimality under correlation ambiguity in binary decisions: de Oliveira–Ishii–Lin — [arXiv](https://arxiv.org/abs/2106.00088)
- Sharp Fréchet partial-identification bounds for classifier metrics in weak supervision (no weight optimisation): Maia Polo et al. NeurIPS 2024 — [NeurIPS](https://proceedings.neurips.cc/paper_files/paper/2024/hash/f4c6bec746b0aeca8c2cd15096f1ad1f-Abstract.html)
- Second-order PAC-Bayes bounds over ρ-weighted votes with weight optimisation, but not sharp: Masegosa et al. 2020 — [ar5iv](https://ar5iv.labs.arxiv.org/html/2007.13532); Wu et al. 2021 — [PDF](https://papers.nips.cc/paper_files/paper/2021/file/69386f6bb1dfed68692a24c8686939b9-Paper.pdf); Germain et al. 2015, with sharpness only in the continuous relaxation — [ar5iv](https://ar5iv.arxiv.org/html/1503.08329)

**What reviewers will expect (from what comparable papers use)**
- **Classical ensembles and bound comparison:** UCI/LibSVM datasets with random forests and heterogeneous ensembles, as in Wu et al. 2021 — [PDF](https://papers.nips.cc/paper_files/paper/2021/file/69386f6bb1dfed68692a24c8686939b9-Paper.pdf). Bound baselines: FO, TND, CCTND/CCPBB, C-bound, and Theisen/Kim disagreement bounds — [Theisen](https://papers.neurips.cc/paper_files/paper/2023/file/30b6fa308e62ed52180c31ae3ba6bb0a-Paper-Conference.pdf), [Kim](https://arxiv.org/html/2411.00328v1). NN ensembles (CIFAR-type), as in Theisen/Kim.
- **Weak supervision:** WRENCH — [D&B](https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/hash/1c9ac0159c94d8d0cbedc973445af2da-Abstract-round2.html). AwA2 binary tasks, as in Mazzetto. Baselines: Snorkel/MeTaL, FlyingSquid, Dawid–Skene, MV, ALL, Mazzetto's subset selection, B&F HedgeClipper.
- **LLM ensembles and juries:**
  - MMLU/MMLU-Pro, GSM8K/MATH, GPQA, TruthfulQA — [Kim 2026](https://arxiv.org/html/2607.20768v1), [J. Chen 2026](https://arxiv.org/abs/2606.27288).
  - QNLI/QQP binary — [Liu](https://arxiv.org/pdf/2605.03379).
  - Judge benchmarks: UltraFeedback, RewardBench, HelpSteer2, ChaosNLI — [Ai et al.](https://arxiv.org/pdf/2510.01499), [Kohli](https://arxiv.org/html/2605.29800v1), [ROPOLL](https://arxiv.org/pdf/2606.30931).
  - LLM aggregation baselines: uniform MV, accuracy-weighted/Nitzan–Paroush log-odds, best single model, Dawid–Skene, spectral meta-learner, OW/ISP (Ai et al. ICML 2026), Ising aggregation (Balasubramanian et al.), CARE, MI-greedy selection (Turkmen et al.).

### Inferences
**Novelty verdict per idea**
- **(A) is partly novel.** Novel parts:
  - the explicit even/odd (label-flip-symmetric versus antisymmetric) decomposition and the theorem that a self-dual weighted vote's accuracy is an odd-part functional, making unlabeled data provably unable to certify it alone;
  - sharp identified sets for weighted votes from (p_i, a_ij), with closed forms for K = 3 and an LP for general K;
  - minimax deterministic weights, plus comparison to B&F's randomized minimax.

  Not novel:
  - the r_ij identity (Platanios 2014);
  - the existence of sign ambiguity (spectral literature);
  - worst-case unweighted subset MV from (ε_i, d_ij) (Mazzetto 2021);
  - transductive minimax with full unlabeled predictions (B&F 2015).
- **(B) appears novel in ML.** I found no ML paper invoking Boros–Prékopa/Prékopa binomial-moment sharp bounds for MV. The C-bound is sharp only for continuous margins and when C ≤ 2·Gibbs; Wu et al. concede FO can beat second-order bounds; Liu's sharp two-moment bounds assume conditionally i.i.d. calls (de Finetti), which is a strictly smaller feasible set than arbitrary dependence among K fixed classifiers. Main risk: a reviewer may see it as "applying a known OR inequality." Mitigate by:
  - (i) proving the PAC-Bayes bounds are strictly non-sharp, with explicit gaps;
  - (ii) handling weighted votes (non-binomial support);
  - (iii) giving finite-sample, PAC-Bayes-compatible versions;
  - (iv) optimising ρ under the sharp bound and showing better certified accuracy or test error.
- **(C) is largely anticipated as a phenomenon:** Kuncheva 2003 for marginals; Chen et al. 2024 for heterogeneous difficulty; Liu 2026 for oscillating curves; Mazzetto Prop. 2 for MV-of-all worse than the median member; Bay & Yearick for plateaus. A worst-case-under-pairwise-information, heterogeneous-member version with LLM evidence (e.g. adding a same-family model lowers the certified lower bound) is still a reasonable section, but not a headline.

**Recommended experiment plan for reviewers**
1. Synthetic: exact LP-sharp intervals versus FO/TND/CC/C-bound/Theisen/Kim bounds as functions of K and correlation, demonstrating strict non-sharpness.
2. Coverage and validity of certified lower bounds with estimated p_i (small labeled n) and a_ij (large unlabeled m), with finite-sample confidence sets.
3. Minimax weights versus uniform, log-odds, PAC-Bayes-optimised, B&F and model-based (DS/Snorkel/Ising/CARE) weights, reporting both worst-case-certified accuracy and realised test accuracy.
4. LLM jury and multi-model voting on judge benchmarks (RewardBench/UltraFeedback/ChaosNLI) and QA (MMLU-Pro/GSM8K binary-ised or correctness-coded), including same-family redundant models to show minimax weights down-weighting redundancy.
5. An add-or-drop-a-model decision rule with a certified guarantee.

A required ablation is a misspecified-dependence stress test, where CI and Ising methods break but the robust weights do not.

**Framing risks to pre-empt**
- (a) Binary only. LLM tasks are multi-class, so the label-flip symmetry argument changes; either restrict to correctness-coded binary decisions or extend.
- (b) Self-duality requires symmetric weights and no abstention.
- (c) The unlabeled agreement P(f_i = f_j) equals P(Z_i = Z_j) only for binary labels. For multiclass, "both wrong but different answers" breaks the identity; Kim et al. ICML 2025 find 60% same-wrong-answer, so the gap is material.

### Gaps
- I did not verify whether any NeurIPS 2025 or ICLR 2026 accepted paper (as opposed to arXiv) covers sharp weighted-vote bounds. OpenReview was not searched directly. A direct OpenReview search for "majority vote" with "sharp" or "partial identification" in the 2025–26 cycles is recommended before submission.
- Unverified venues: Lacasse et al. 2007; Parisi et al. PNAS 2014; Arachie & Huang AAAI 2019; Verga et al. 2024; Guo et al.; Huang–Xia–Kong; Denisov-Blanch et al.; and all 2026 arXiv items.
- Not read in full, but potentially overlapping: Corrada-Emmanuel NTQR papers (arXiv 2312.05392, 2509.08593, 2510.00821); "On the generalization of the C-bound" (arXiv 1408.1336); Biggs, Zantedeschi, Guedj, "On margins and generalisation for voting classifiers" (NeurIPS 2022).
