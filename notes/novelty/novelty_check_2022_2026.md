# Novelty check (2022–2026): weighted voting / ensemble aggregation when marginal accuracies are known but dependence is uncertain

Scope and conventions:
- This note checks the candidate topics extracted from the four earlier notes in this folder:
  - `ml_bounds_unlabeled_llm.md`
  - `probability_or_robust_aggregation.md`
  - `social_choice_and_fusion_weights.md`
  - `unlabeled_agreement_sharp_bounds.md`
- It does not re-derive anything.
- Verification tags:
  - **[FT]**: checked against full text or a full HTML/PDF rendering, this session or the earlier notes.
  - **[ABS]**: abstract or landing page only.
  - **[PRIOR-NOTE]**: taken from the earlier notes in this folder (which give their own sources) and not re-fetched this session.
- Search date: 13 Oct 2026. "Not found" means absent from about 45 targeted searches and fetches across arXiv, PMLR, NeurIPS, JMLR, OpenReview listings, IDEAS and Springer. It is not proof of absence.

**Candidate claims extracted from the notes:**
- (a) **Max–min weighted votes.**
  - a1: with marginals only, the max–min rule is the best single classifier (a "dictator"). The notes derive this from a covering-number formula.
  - a2: with marginals plus pairwise information (equivalently, labeled accuracies plus unlabeled pairwise agreement), K=3 has a closed form. Use majority iff Σ_{i<j}T_ij < e_min, where T_ij = P(both i and j wrong) and e_min is the smallest individual error; otherwise use the best single classifier. Randomising never helps.
  - a3: K=5 phase diagrams. In the notes' numerical experiments, the max–min rule was always either a single classifier or a majority of 3 on positively dependent instances. Signed (negative-weight) rules were sometimes optimal for K=3.
  - a4: general-K complexity.
- (b) **Partial identification.**
  - b1: an even/odd identification theorem. Unlabeled data identify the label-flip-symmetric ("even") part of π, the joint law of the correctness pattern; self-dual vote accuracy depends only on the odd part. The identified set is E(m,p), and only single-classifier accuracies are point-identified.
  - b2: sharp "pair-splitting" LP intervals for any weighted vote, plus a K=3 closed form.
  - b3: a sharp replacement for PAC-Bayes second-order bounds. As K→∞ it is min{C-bound, 3μ−2s}, i.e. L(MV) ≤ E_ρ[L]+E_{ρ²}[D] (individual error plus pairwise disagreement). For finite K it has a Boros–Prékopa-type closed form, and for weighted votes it is a moment LP over the possible weighted error sums V_ρ.
  - b4: max–min and minimax-regret weights built from these sets.
- (c) **Applications.** Correlated LLM ensembles, LLM-judge panels, and weak-supervision label models when the assumed dependence structure is wrong. A further option is non-monotonicity in K (adding voters can hurt) as a showcase.

---

## Q1. Candidate (a): max–min optimal weighted-vote rules over Fréchet-type ambiguity sets (K=3, K=5, general K)

### Takeaway
- **a1 (marginals only) is already known.** Do not claim it.
- **a2 (K=3 with pairwise information) is essentially contained in Mazzetto et al. (AISTATS 2021).** They give a closed form for the worst-case error of a 3-voter majority from error rates and pairwise disagreements, and they choose among odd-sized unweighted sub-majorities by worst case. For K=3 every self-dual monotone rule is a single classifier or the majority, so a2 follows in a few lines.
- **a3 is open but at high risk of being "Mazzetto with weights."** The notes' own experiments mostly selected rules (single classifier or 3-majority) that are already in Mazzetto's candidate class. A K=5 phase diagram only adds novelty in the regimes where a genuinely weighted rule, such as (2,1,1,1,0) or (3,1,1,1,1), or a signed rule wins.
- **a4 (complexity) is partly known.** Evaluating the worst case of one rule is NP-hard by reductions that already exist. Hardness or approximability of the *outer* max over weights under pairwise information is not addressed anywhere I found.
- de Oliveira–Ishii–Lin explicitly leave "known pairwise correlations" as an open extension. That is the cleanest citation for why the pairwise max–min problem is open.

### Cited Findings
**Marginals only: the best single source is max–min optimal (kills a1 as a contribution)**
- de Oliveira, Ishii & Lin, "Robust Aggregation of Correlated Information," arXiv 2106.00088 (latest version dated Sept 2026; earlier "Robust Merging of Information," ACM EC 2021, DOI 10.1145/3465456.3467533).
  - Theorem 1: with two states and two actions, a best-source strategy is robustly optimal against all joint laws with the given marginals. [FT, prior note] — [arXiv](https://arxiv.org/html/2106.00088)
  - Section 7.1.2 "Partial Knowledge of Correlations": known block-joint laws reduce to the marginal case. For pairwise knowledge (e.g. P1–P2 and P2–P3 known) the authors write that "our analysis does not immediately extend to other, more complex situations" and leave general restrictions to future research. [FT this session] — [arXiv html](https://arxiv.org/html/2106.00088)
- Arieli, Babichenko, Talgam-Cohen & Zabarnyi, "A random dictator is all you need," *AEJ: Microeconomics* 17(1):66–96 (2025); arXiv 2302.03667. With symmetric agents, marginals known and correlation adversarial, a random dictator is uniquely optimal for many agents under maximin, regret and approximation ratio; simple majority is "far from optimal." [ABS + prior note] — [IDEAS](https://ideas.repec.org/a/aea/aejmic/v17y2025i1p66-96.html)
- Mazzetto, Sam, Park, Upfal & Bach, "Semi-supervised aggregation of dependent weak supervision sources with performance guarantees," AISTATS 2021, PMLR 130:3196–3204. Prop. 3: with only labelers' outputs and error rates, max over consistent laws of min over all functions f of ε(f∘ℓ_S) = min{ε_1,…,ε_n}, so no aggregator beats the best labeler. [FT this session] — [NSF PAR](https://par.nsf.gov//servlets/purl/10273522); [PMLR](https://proceedings.mlr.press/v130/mazzetto21a.html)
- Balsubramani & Freund, COLT 2015 (PMLR 40). Without the unlabeled prediction matrix F, the minimax solution reduces to the classifier with the best error guarantee. [FT, prior note] — [arXiv](https://arxiv.org/pdf/1503.01811)

**Marginals plus pairwise information: closest prior art to a2 and a3**
- Mazzetto et al. 2021 [FT this session]:
  - **Prop. 5** gives the worst-case majority error of 3 labelers {i,j,k} in closed form from ε and pairwise disagreements D:
    - q = ½·max{ε_i+ε_j−D_ik−D_jk, ε_j+ε_k−D_ij−D_ik, ε_i+ε_k−D_ij−D_jk, 0};
    - worst-case error = ε_i+ε_j+ε_k − (D_ij+D_jk+D_ik)/2 − 2q;
    - in the equal-error case, 3(ε−d/2) − 2max{ε−d, 0}.
  - The worst-case error of a general subset is an LP over correctness patterns {0,1}^k, with constraints: individual errors, pairwise disagreements, normalisation and nonnegativity.
  - The algorithm picks the odd subset I minimising this worst case. For large n it uses the size-3 closed form ("PGMV") or greedy growth by pairs ("PGMV-P/-D"). Singletons are allowed.
  - Candidates are **unweighted** majorities only. The paper claims no minimax optimality beyond this class and gives no finite-sample confidence sets.
  - Experiments: 45 binary AwA2 tasks.
  - Sources: [NSF PAR](https://par.nsf.gov//servlets/purl/10273522)
- Kaniovski, *Theory and Decision* 69(3):453–468 (2010), and Zaigraev & Kaniovski, *Theory and Decision* 72(1):89–112 (2012). These give simple-majority bounds under exchangeable correlation, with homogeneous min M = (n(2p−1)+1)/(n+1). They do not optimise weights. [PRIOR-NOTE] — [preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_minmax.pdf)
- Althöfer & Thiele, "A Condorcet jury theorem for couples," *Theory and Decision* 81(1):1–15 (2016), DOI 10.1007/s11238-015-9521-0. Tight lower and upper bounds for majority over N independent couples, where the two agents in each couple are correlated. The lower bound equals the majority of N independent agents. This is a structured-dependence result, not max–min weighting. [ABS] — [IDEAS](https://ideas.repec.org/a/kap/theord/v81y2016i1d10.1007_s11238-015-9521-0.html)
- Pan, Chen & Kong, "Robust decision aggregation with second-order information," WWW 2024 (arXiv 2311.14094). Two experts only. Over all information structures, second-order information gives no benefit (regret 0.5, "follow the first expert"). [PRIOR-NOTE] — [arXiv](https://arxiv.org/html/2311.14094v1)
- Guo & Kong, "Robust Decision Aggregation with Adversarial Experts," arXiv 2403.08222 (v2 6 Feb 2025). Binary states with marginally symmetric experts, some of them adversarial. A (randomised) truncated mean is optimal. The ambiguity concerns adversarial reports, not correlation among honest experts with known accuracies. [ABS] — [arXiv](https://arxiv.org/abs/2403.08222)
- Chen, Peng & Tang, "Prior-agnostic robust forecast aggregation," arXiv 2604.24517 (2026). Over general correlated structures the minimax regret is 1/4 and is trivial; positive results need conditional independence (CI) or Blackwell order. [PRIOR-NOTE] — [arXiv](https://arxiv.org/pdf/2604.24517)

**Complexity (a4)**
- Kreinovich & Ferson, *IJAR* 41:331–342 (2006). Best-possible bounds on P(Σx_i<0) for two-point marginals are NP-hard, by reduction from Partition. [PRIOR-NOTE] — [tech report](https://www.cs.utep.edu/vladik/2004/tr04-25a.pdf)
- Padmanabhan, Ahipasaoglu, Ramachandra & Natarajan, "Extremal probability bounds in combinatorial optimization," arXiv 2109.01591 (SIAM, DOI 10.1137/21M1442504).
  - Tight upper bounds given marginals: weakly NP-hard, with a pseudo-polynomial algorithm.
  - Lower bounds: strongly NP-hard even for Bernoulli variables.
  - [PRIOR-NOTE] — [arXiv](https://arxiv.org/pdf/2109.01591)
- Pitowsky, *Math. Programming* 50:395–414 (1991). Membership in the correlation polytope (are given p_i, p_ij consistent?) is NP-complete. [PRIOR-NOTE]
- Kaski, Mannila & Mohapatra, "Optimal union probability interval is NP-hard," ESA 2026, LIPIcs 388, 129. [PRIOR-NOTE] — [Dagstuhl](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ESA.2026.129)
- Tang & Zhang, "Algorithmic Expert Aggregation," arXiv 2607.08744 (2026). Not robust-over-correlation. Its deterministic aggregation version is NP-hard even with two experts, while the randomised version is polynomial. This is a useful analogy for a "deterministic vs randomised" gap. [FT this session] — [arXiv](https://arxiv.org/html/2607.08744v1)
- Neoh, Teh, Chooi, Goldberg & Tambe, "Efficient ensemble selection from binary and pairwise feedback," arXiv 2605.09588 (2026). Coverage-based committee selection with (1−1/e) guarantees and a PTAS-but-no-EPTAS result. It has no majority-vote accuracy and no weights. [FT this session] — [arXiv](https://arxiv.org/html/2605.09588)

### Inferences
- **a2 novelty risk: HIGH.**
  - For K=3 the only self-dual monotone rules are the three single classifiers and the majority. Mazzetto Prop. 5 gives the majority's worst case, and a single classifier's accuracy is identified, so "max–min = majority iff Σ T_ij < e_min" is a short corollary.
  - Remaining contributions: the T-form; the symmetric phase boundary ρ*(p); and the observation that signed rules can win once non-monotone rules are allowed. Signed rules fall outside Mazzetto's class and need exact accuracies, not lower bounds (contrast Balsubramani–Freund's σ ≥ 0).
- **a3 novelty risk: MEDIUM–HIGH.**
  - New: weighted classes and signed rules, and a characterisation theorem ("when is the max–min rule a committee of size ≤ 3?").
  - Risk: if the main finding is "the max–min rule is a single classifier or a 3-majority," reviewers will say Mazzetto's PGMV already searches exactly that class.
  - What would separate the work is a *proved* structural theorem, e.g. "under the pairwise set the max–min rule is always attained by a committee of size ≤ 3", or a counterexample family. Such a theorem would also explain *why* Mazzetto's size-3 heuristic works. That is a genuinely interesting claim. I found no source proving or refuting it.
- **a4 novelty risk: MEDIUM.**
  - NP-hardness of evaluating one rule is inherited from Kreinovich–Ferson, Padmanabhan et al. and Pitowsky.
  - Not found anywhere:
    - hardness or approximation of the outer max–min over weights;
    - fixed-parameter or tree-structured tractability for weighted thresholds with pairwise data (extending Natarajan–Padmanabhan's tree LP from counts to weighted sums);
    - a dictator-to-committee phase-transition theorem as the ambiguity radius shrinks.
- The phrase "distributionally robust ensemble weighting" with a Fréchet or pairwise ambiguity set is essentially unused in ML venues. That makes positioning easy, but reviewers from the DRO community will expect a link to Wasserstein or φ-divergence sets.

### Gaps
- I did not read Narasimhamurthy (TPAMI 2005), the LP bounds for majority with pairwise information. It may already contain the K=3 worst-case closed form, pre-dating Mazzetto.
- I did not read Berg (1994), "Evaluation of some weighted majority decision rules under dependent voting" (*Math. Social Sciences* 28(2):71–83), which the prior note flagged. It must be read before claiming novelty for weighted rules under dependence.
- I found no 2022–2026 follow-up to Mazzetto et al. 2021 that adds weights. Searches for "Mazzetto … weighted … follow-up" returned only the original.

---

## Q2. Candidate (b): sharp partial-identification bounds from labeled marginals plus unlabeled agreement, robust weights, and comparison with PAC-Bayes C-bound, tandem and disagreement bounds

### Takeaway
- **b1 (even/odd identification)** is not stated in this form anywhere I found. The ingredients are known: the Platanios agreement identity, the sign ambiguity noted in the spectral literature, and Balsubramani–Freund's game. A clean theorem is publishable as part of a paper, but not as its headline. Risk is MEDIUM.
- **b2 (sharp intervals for a given weighted vote from (p, pairwise agreement) or (p, unlabeled law))** is a direct weighted generalisation of Mazzetto's LP. The structural "pair-splitting LP" and its dual certificate appear new. Risk is MEDIUM.
- **b3 (sharp replacement for the PAC-Bayes bounds)**: much of it is already in print.
  - **The K→∞ formula min{C-bound, 3μ−2s} is exactly Theorem 4.3 of Liu (arXiv 2605.03379, 2026).** There it is stated for the mixing law of the per-example correctness probability Q, which is mathematically the same moment problem on [0,1]. It is also a classical Markov–Krein/Karlin–Studden moment bound.
  - The finite-K uniform case is almost certainly a special case of Boros–Prékopa (1989, MOR 14(2):317) two-binomial-moment bounds.
  - What remains new: the *translation* into the PAC-Bayes language, L(MV_ρ) ≤ E_ρ[L] + E_{ρ²}[D]; the demonstration that the C-bound, tandem and first-order bounds are never sharp for finite K; the weighted (support V_ρ) version; and per-classifier (p, r) bounds.
  - I did not find "E[L]+E[D]" in Lacasse et al. 2006, Germain et al. 2015, Laviolette et al. 2014/15, Wu et al. 2021, Theisen et al. 2023 or Feofanov et al. 2024. Theisen et al. do have the mirror-image *lower* bound L(MV) ≥ E[L] − E[D].
  - Overall risk is HIGH if presented as a new inequality and MEDIUM if presented as "sharpness theory for majority-vote certificates."
- **b4 (robust weights)** inherits the issues of Q1.

### Cited Findings
**Identification ingredients (relevant to b1)**
- Platanios, Blum & Mitchell, UAI 2014. Pairwise agreement identity a_ij = 1 − e_i − e_j + 2e_ij; the system is underdetermined without independence. [FT, prior note] — [PDF](https://www.auai.org/uai2014/proceedings/individuals/313.pdf)
- Parisi et al. (PNAS 2014) and Jaffe, Nadler & Kluger (AISTATS 2015). Spectral identification under CI, up to a ±1 sign ambiguity. [FT, prior note] — [ar5iv](https://ar5iv.arxiv.org/html/1407.7644)
- Corrada-Emmanuel, "Logical Consistency Between Disagreeing Experts and Its Role in AI Safety," arXiv 2510.00821 (1 Oct 2025). [FT this session] — [arXiv html](https://arxiv.org/html/2510.00821v1)
  - Computes the set of evaluations logically consistent with observed agreements and disagreements, via an LP over integer response counts, with "axioms" (linear equalities) obtained from label symmetry.
  - Works mainly with M=1 per-classifier counts; Appendix C adds pairwise (M=2) counts.
  - States that "M>2 summaries must lie within the product of M=1 sets" but does not work out triples or general K.
  - Gives no bound on ensemble or majority-vote accuracy and no even/odd decomposition.
- Corrada-Emmanuel, "No-Knowledge Alarms for Misaligned LLMs-as-Judges," arXiv 2509.08593 (10 Sept 2025). LP over integer counts; alarms that flag, with no false positives, that some judge violates a user-set grading floor. [ABS] — [arXiv](https://arxiv.org/abs/2509.08593)
- Corrada-Emmanuel, "Algebraic Evaluation Theorems," arXiv 2412.16238. Under error independence, two solutions exist (a label-flip-type ambiguity). [PRIOR-NOTE] — [arXiv](https://arxiv.org/html/2412.16238v2)
- Haber & Wexler, "Unsupervised Estimation of Ensemble Accuracy," arXiv 2311.10940 (preprint "under review"; an mlanthology entry suggests a NeurIPS 2023 workshop). A label-free combinatorial bound on ensemble mistakes; its computation is NP-hard, with a PTAS and a greedy 4/3-approximation. It assumes errors pick the target class uniformly at random, and does not cover weighted votes. [FT this session] — [arXiv](https://arxiv.org/html/2311.10940v2); [mlanthology](https://mlanthology.org/neuripsw/2023/haber2023neuripsw-unsupervised)

**Sharp evaluation of a vote's accuracy (relevant to b2)**
- Mazzetto et al. AISTATS 2021. LP worst case for unweighted subset majority from (ε, D), or from the full output distribution Q. [FT this session] — [NSF PAR](https://par.nsf.gov//servlets/purl/10273522)
- Maia Polo, Maity, Yurochkin, Banerjee & Sun, "Weak supervision performance evaluation via partial identification," NeurIPS 2024 (arXiv 2312.04601). [FT this session; venue from the prior note] — [arXiv](https://arxiv.org/pdf/2312.04601); [NeurIPS](https://proceedings.neurips.cc/paper_files/paper/2024/hash/f4c6bec746b0aeca8c2cd15096f1ad1f-Abstract.html)
  - Fréchet bounds on a classifier's accuracy or F1, using P(X,Z) from unlabeled data and P(Y|Z) from a label model.
  - No weighted votes or weight optimisation; dependence among labeling functions enters only through the label model.
- Feofanov, Devijver & Amini, "Multi-class probabilistic bounds for majority vote classifiers with partially labeled data," *JMLR* 25 (2024), 1–47. [FT this session] — [JMLR](https://jmlr.org/papers/volume25/23-0121/23-0121.pdf)
  - A transductive LP-based majority-vote error bound from votes and the Gibbs conditional risk, plus C-bound variants with imperfect labels.
  - Dependence is handled only through the second margin moment. No weight optimisation, and no sharpness over dependence.
- Shen, Balashankar, Fisch, Madras & Rodrigues, "Robust LLM performance certification via constrained maximum likelihood estimation," arXiv 2604.03257 (2026). Certifies a single LLM's failure rate using one judge, with a small labeled set plus a large judge-labeled set. It assumes i.i.d. data and has no ensemble or voting. [FT this session] — [arXiv](https://arxiv.org/html/2604.03257v1)
- Balsubramani & Freund, COLT 2015 / NeurIPS 2016. A transductive minimax game whose solution is a clipped weighted vote, with no independence assumption; σ ≥ 0 because accuracies enter as lower bounds. [FT, prior note] — [arXiv](https://arxiv.org/pdf/1503.01811)

**PAC-Bayes and second-order majority-vote bounds (relevant to b3)**
- Lacasse, Laviolette, Marchand, Germain & Usunier, "PAC-Bayes bounds for the risk of the majority vote and the variance of the Gibbs classifier" (NIPS 2006 per the PDF file name `nips06a`; venue not printed in the text). [FT this session] — [PDF](http://www.pascalgermain.info/publis/nips06a.pdf)
  - Bounds: R(B_Q) ≤ 2R(G_Q), and the Cantelli bound C_Q = Var(W_Q)/(Var(W_Q) + (½−R(G_Q))²), where B_Q is the majority vote and G_Q the Gibbs classifier.
  - Calls Cantelli "the tightest possible upper bound" given mean and variance, which holds for general real random variables, not for [0,1]-valued or lattice ones.
  - Prop. 7: degradation as d_Q→½ is unavoidable for bounds that use only R(G_Q) and Var(W_Q).
  - No bound of the form R(G) + d.
- Germain, Lacasse, Laviolette, Marchand & Roy, *JMLR* 16:787–860 (2015). C-bound; Prop. 13: sharp iff μ₂ ≤ μ₁; MinCq. [FT, prior note] — [ar5iv](https://ar5iv.arxiv.org/html/1503.08329)
- Laviolette, Morvant, Ralaivola & Roy, "On the generalization of the C-bound to structured output ensemble methods," arXiv 1408.1336 (technical report v2, 2015). Multiclass and multilabel C-bounds; no tightness theorem, no finite-voter version, no R(G)+d form, no moment-problem analysis. [FT this session] — [alphaxiv](https://www.alphaxiv.org/abs/1408.1336.md)
- Masegosa, Lorenzen, Igel & Seldin, NeurIPS 2020 (tandem bound), and Wu, Masegosa, Lorenzen, Igel & Seldin, NeurIPS 2021 (Chebyshev–Cantelli / PAC-Bayes-Bennett). Wu et al. concede that "in the high-error-correlation regime the first order oracle bound is tighter than the second order oracle bounds." [FT, prior note] — [NeurIPS 2021 PDF](https://papers.nips.cc/paper_files/paper/2021/file/69386f6bb1dfed68692a24c8686939b9-Paper.pdf)
- Wu & Seldin, "Split-kl and PAC-Bayes-split-kl inequalities for ternary random variables," NeurIPS 2022 (arXiv 2206.00706). A concentration tool for tandem/disagreement estimates, not a new oracle bound. [ABS] — [arXiv](https://arxiv.org/abs/2206.00706)
- Theisen, Kim, Yang, Hodgkinson & Mahoney, "When are ensembles really effective?", NeurIPS 2023 (arXiv 2305.12313). [FT this session (supplement)] — [supplement](https://proceedings.nips.cc/paper_files/paper/2023/file/30b6fa308e62ed52180c31ae3ba6bb0a-Supplemental-Conference.pdf)
  - Thm 3 gives a lower bound without any competence assumption: **L(MV) ≥ E_ρ[L] − E_{ρ²}[D]**.
  - Upper bound under competence: L(MV) ≤ 4(K−1)/K·(E[L] − ½E[D]), where K here is the number of classes; for binary tasks this is 2E[L] − E[D].
  - Tightness examples (App. C) exist only for the 2E[L] bound.
- Kim, Hodgkinson, Theisen & Mahoney, "How many classifiers do we need?", NeurIPS 2024 (arXiv 2411.00328). η-polarization bounds with an N/(N+1) finite-N factor; no sharpness. [FT, prior note] — [arXiv](https://arxiv.org/html/2411.00328v1)
- Hauptvogel & Igel, arXiv 2406.05469 (2024/25): tandem-bound weighting of deep ensembles. Bastian, Leblanc, Germain, Habrard, Metzler, Morvant & Viallard, "On the disintegration of the stochastic majority vote," arXiv 2609.16803 (15 Sept 2026): deterministic-vote PAC-Bayes bounds. Neither claims sharpness or uses disagreement or unlabeled data in the abstract. [ABS] — [arXiv](https://arxiv.org/abs/2609.16803)
- **Liu (York University), "Two Calls, Two Moments, and the Vote-Accuracy Curve of Repeated LLM Inference," arXiv 2605.03379 (2026).** [FT this session] — [arXiv PDF](https://arxiv.org/pdf/2605.03379)
  - **Thm 4.3 (infinite votes)** gives sharp bounds on the infinite-vote accuracy V_∞ given μ = E[Q] and ν = E[Q²], where Q is the per-example probability that a call is correct. U_∞ = G(μ,ν), with:
    - G = (ν−μ²)/(ν+¼−μ) if μ ≤ ½ and ν ≤ μ/2;
    - G = min{1, 3μ−2ν} otherwise.
    - L_∞ = 1 − G(1−μ, 1−2μ+ν).
  - Translated to error (W = 1−Q), this is **identical in form** to the notes' "Theorem 5" (min{C-bound, 3μ−2s}).
  - Thm 4.2 gives three-vote closed forms L_1 and U_1. Calls are assumed conditionally i.i.d. given Q.
  - Cites Germain et al. 2015 for the C-bound and PAC-Bayes, but not Masegosa, Prékopa or Kuncheva.
  - Weighted votes and heterogeneous models are out of scope (heterogeneity appears only empirically, Appendix F.1).
- Boros & Prékopa, "Closed form two-sided bounds for probabilities that at least r and exactly r out of n events occur," *MOR* 14(2):317–342 (1989). Subasi, Subasi & Prékopa (RUTCOR RRR 4-2006) confirm that Boros–Prékopa give closed forms for the at-least-r probability "when m ≤ 4" moments. [ABS/secondary] — [RRR 4-2006](https://rutcor.rutgers.edu/Prekopa/pdf/4_2006.pdf); [INFORMS](https://pubsonline.informs.org/doi/fpi/10.1287/moor.14.2.317)

### Inferences
- **b3 has a concrete novelty problem.**
  - The headline "sharp bound min{C, μ+E[D]}" is the same two-moment problem on [0,1] that Liu solved in 2026, and it is a classical Markov–Krein extremal problem.
  - The PAC-Bayes community apparently never wrote the linear branch μ + E[D] (verified absent from the six papers listed above). It is a nice *remark*, the mirror of Theisen's lower bound E[L]−E[D], but it will be read as "known moment theory applied."
  - Defensible novelty:
    - (i) for a *fixed finite* weighted ensemble the C-bound is sharp on 0% of the grid, so every published second-order bound is strictly non-sharp. This needs a theorem, not just a grid.
    - (ii) the weighted support V_ρ reduction and its closed form;
    - (iii) the per-classifier (p, r) and (p, unlabeled-law) bounds, which beat any averaged-moment bound by a large margin (the notes' K=5 log-odds example: 0.313 vs 0.187);
    - (iv) PAC-Bayes-kl finite-sample versions plus optimisation over ρ.
- **b1 + b2 together are the most original theory in the folder.** The cleanest statement is: "the identified set is E(m,p); only single-classifier accuracies are point-identified; the worst-case adversary orients all but ≤ K antipodal pairs." I found no paper with this statement. The Corrada-Emmanuel line is the most likely overlap and stops at pairs with no ensemble-vote bounds.
- The finite-sample message from the notes' experiments (the dominant uncertainty is labeled p̂, not identification) is not in Mazzetto, Maia Polo or Feofanov. It is a useful practical contribution.

### Gaps
- I did not obtain the full text of Boros–Prékopa 1989 to confirm that the notes' "Theorem 6" coincides with their m=2 formula.
- I did not search exhaustively for whether any 2025–26 NeurIPS/ICLR/AISTATS paper contains the bound "R(MV) ≤ R(G) + d". An OpenReview full-text search is recommended. Six core papers were checked and none contain it.
- arXiv 2505.21972 appeared in a search on partial identification for LLM judges but could not be fetched (robots-disallowed); its content is unknown.

---

## Q3. Candidate (c): LLM-voting and LLM-judge ensembles, and weak supervision under misspecified dependence

### Takeaway
- **The empirical motivation is very strong and recent.** Several 2025–26 papers document correlated LLM and judge errors:
  - Kim et al. ICML 2025; Goel et al. ICML 2025.
  - Kohli 2026: 9 judges ≈ 2 effective votes.
  - J. Chen 2026: a co-failure ceiling; pairwise correlation cannot identify the all-wrong rate.
  - Kim 2026: majority voting beats the best member in only about 10% of size-3 LLM subsets.
- Every aggregation *method* with theory assumes a dependence model:
  - CI: Ai et al. ICML 2026.
  - Ising model: Balasubramanian et al. 2026.
  - Latent confounders: CARE, ICML 2026.
  - Gaussian copula: Turkmen et al. 2026.
  - Bradley–Terry with product likelihood: Qian et al. 2026.
  - Equicorrelation / ε-contamination: ROPOLL 2026.
- Repeated-sampling theory assumes exchangeability (Liu 2026 ×2; Chen et al. NeurIPS 2024).
- **I found no paper giving assumption-free (Fréchet/pairwise-ambiguity) certified intervals or max–min weights for heterogeneous LLM or judge panels.** As an *application track* for (a)+(b) this is novel. As a stand-alone "LLMs are correlated" paper it is not.
- For weak supervision, misspecified dependence has been studied, but no paper gives sharp worst-case weighted-vote guarantees:
  - Cachay et al. 2021: structure over-specification bounds;
  - Mazzetto et al. 2021: worst-case unweighted subsets;
  - Maia Polo et al. NeurIPS 2024: Fréchet evaluation bounds.

### Cited Findings
**Correlated errors (empirical)**
- Kim, Garg, Peng & Garg, "Correlated errors in large language models," ICML 2025, PMLR 267. When two models both err, they agree on the same wrong answer 60% of the time on one leaderboard dataset. [PRIOR-NOTE] — [PMLR](https://proceedings.mlr.press/v267/kim25e.html)
- Goel et al., "Great models think alike and this undermines AI oversight," ICML 2025, PMLR 267. [PRIOR-NOTE] — [PMLR](https://proceedings.mlr.press/v267/goel25b.html)
- Kohli, "Nine judges, two effective votes," arXiv 2605.29800 (May 2026). [FT this session] — [arXiv](https://arxiv.org/html/2605.29800v1)
  - A "Condorcet null" simulation under CI (per-judge, per-difficulty-bin confusion matrices); effective votes via Kish n_eff = k/(1+(k−1)φ̄) and k/λ_max.
  - Oracle Dawid–Skene and accuracy-weighted voting close at most 11% of the Condorcet gap.
  - **No worst-case bounds or robust weighting.** The paper recommends diversifying the panel instead.
- J. Chen, "When does combining language models help? A co-failure ceiling on routing, voting, and mixture-of-agents across 67 frontier models," arXiv 2606.27288 (25 June 2026). [ABS-level page this session] — [arXiv](https://arxiv.org/abs/2606.27288)
  - Accuracy ≤ 1 − β for member-selecting policies, where β is the all-wrong rate.
  - Equal marginals and pairwise correlations can give different β; a Gaussian copula underprices β about 2.5×. β is certified with Clopper–Pearson.
  - No sharp or LP bounds for votes and no weight selection.
- Kim (Aidentyx), "Are diversity metrics measuring diversity?", arXiv 2607.20768 (2026). [PRIOR-NOTE] — [arXiv](https://arxiv.org/pdf/2607.20768)
- Ali, "Quantifying Diversity of Thought: A Predictive Law of Weighted LLM Ensemble Lift," arXiv 2607.17384 (mid-2026). [FT this session] — [HF papers](https://huggingface.co/papers/2607.17384.md)
  - A two-model exact lift decomposition with a correlation term φ_adj.
  - An oracle-selection ceiling min(q, 1−p).
  - Three-model results only empirical; no general worst-case bounds.
- Ding, "When LLMs agree, are they right?", arXiv 2607.08065 (v2 28 July 2026). Purely empirical; no aggregation theory. [FT this session] — [arXiv](https://arxiv.org/html/2607.08065)
- Sasahara, Naito & Fujie, "A latent dimension of Condorcet's jury theorem for multiple AI advisers," arXiv 2609.14438 (13 Sept 2026). [FT this session] — [arXiv](https://arxiv.org/html/2609.14438v1)
  - Main text assumes independent advisers of equal accuracy; dependence appears only in Gaussian-copula simulations.
  - A "dissent deficit" diagnostic. No worst-case bounds or weights.

**LLM aggregation methods with theory (all model-based)**
- Ai, Pan, Simchi-Levi, Tambe & Xu, "Beyond majority voting: LLM aggregation by leveraging higher-order information," ICML 2026, PMLR 306. CI given the label; no minimax guarantee. [PRIOR-NOTE] — [ICML](https://icml.cc/virtual/2026/poster/63160)
- Zhao, Shin, Huang, Namburi & Sala, "CARE: Confounder-aware aggregation for reliable LLM evaluation," ICML 2026, PMLR 306:162114–162150. Latent-confounder model; identifiability and finite-sample recovery under that model; no worst-case guarantees. [ABS this session] — [PMLR](https://proceedings.mlr.press/v306/zhao26aq.html)
- Balasubramanian, Podkopaev & Kasiviswanathan, "Dependence-aware label aggregation for LLM-as-a-judge via Ising models," arXiv 2601.22336 (2026). [PRIOR-NOTE] — [arXiv](https://arxiv.org/pdf/2601.22336)
- Turkmen, Buyukates & Bastopcu, "Don't always pick the highest-performing model," arXiv 2602.08003 (2026). Gaussian copula. [PRIOR-NOTE] — [arXiv](https://arxiv.org/pdf/2602.08003)
- Qian, Sun, Gales & Knill, "Who can we trust? LLM-as-a-jury for comparative assessment," arXiv 2602.16610 (18 Feb 2026). [FT this session] — [arXiv](https://arxiv.org/html/2602.16610v1)
  - BT-σ: per-judge discriminators learned without labels.
  - The likelihood is a product over judges, which is implicitly CI (my reading of Eq. 11; the authors do not state it).
  - No robustness guarantees.
- Acharya, Pan & Verkhovsky (AWS), "ROPOLL," arXiv 2606.30931 (2026). [PRIOR-NOTE] — [arXiv](https://arxiv.org/pdf/2606.30931)
- Exchangeable repeated-sampling theory:
  - Liu, arXiv 2605.03379 and 2605.05592 (2026) [FT 2605.03379 this session];
  - Chen et al., "Are more LLM calls all you need?", NeurIPS 2024 [PRIOR-NOTE];
  - Mattei & Garreau, "Are ensembles getting better all the time?", *JMLR* 26 (2025). Exchangeable predictions, unweighted averaging; ensembles can eventually worsen under non-convex losses or classification error when the asymptotic prediction is wrong. [FT this session] — [JMLR](https://jmlr.org/papers/volume26/24-0408/24-0408.pdf)

**Weak supervision under misspecified dependence**
- Cachay, Boecking & Dubrawski, "Dependency structure misspecification in multi-source weak supervision models," arXiv 2106.10302 (2021; WeaSuL workshop proceedings per arXiv 2107.03690). Structure over-specification can hurt; theoretical bounds on that modelling error. [ABS] — [arXiv](https://arxiv.org/abs/2106.10302v1)
- The Snorkel line: Ratner 2016; Bach ICML 2017; Varma ICML 2019; Fu et al. ICML 2020 (FlyingSquid). Consistency under a correct or sparse graph; no worst-case guarantees. [PRIOR-NOTE] — [FlyingSquid ar5iv](https://ar5iv.labs.arxiv.org/html/2002.11955)
- Arachie & Huang, *JMLR* 22 (2021); Mazzetto et al. ICML 2021. Adversarial label learning with per-source constraints only. [PRIOR-NOTE] — [JMLR](https://www.jmlr.org/beta/papers/v22/20-537.html)

### Inferences
- **An LLM-judge / LLM-voting application of (a)+(b) has LOW–MEDIUM novelty risk** as long as the theory is the contribution and the LLM experiments are the testbed. The 2026 literature has identified the problem (Kohli; J. Chen; Kim) but answers it with either a dependence *model* or an *advice to diversify*. "Certified, assumption-free intervals plus max–min weights from a small labeled set plus agreement on a large unlabeled pool" is unoccupied.
- **Technical caveat that reviewers will raise.** LLM tasks are multiclass, and binary agreement P(f_i = f_j) differs from P(Z_i = Z_j) when two models are both wrong but give different answers. Kim et al. report 60% same-wrong-answer, so the gap is material. The theory must either (i) restrict to binary or pairwise-preference judge tasks (RewardBench, UltraFeedback pairs, ChaosNLI-binarised), or (ii) extend the identification to multiclass, where unlabeled data reveal more than the even part. (ii) is itself a novel and harder theorem.
- **"Non-monotone in K" as a headline: HIGH risk.** It is covered by Chen et al. 2024, Liu 2026 ×2, Mattei & Garreau 2025, Mazzetto Prop. 2, Kuncheva 2003 and Bax 2021. Keep it as a corollary under pairwise information with heterogeneous members.
- **A pure weak-supervision framing: MEDIUM risk.** Maia Polo (NeurIPS 2024) owns "partial identification for WS evaluation," and Mazzetto 2021 owns "worst-case subset majority." The new piece would be weighted votes plus abstention handling (WRENCH labeling functions abstain), which none of these handle under Fréchet ambiguity.

### Gaps
- Venue status of most 2026 arXiv items is unknown: Kohli, J. Chen, Liu, Ali, Sasahara, Qian, Balasubramanian, Turkmen, ROPOLL, Shen.
- I did not check NeurIPS 2026 accepted papers. Notifications were around late Sept 2026, and the proceedings were not indexed in my searches.

---

## Q4. Social-choice side: robust or ambiguity-averse jury theorems and optimal weights under unknown correlation (2015–2026)

### Takeaway
- Social choice has (i) jury theorems under *specific* dependence models (opinion leaders, couples, exchangeable Bahadur), (ii) ambiguity-averse *strategic* voting models, and (iii) "robust voting" papers about preference uncertainty. None of them chooses weights max–min over joint competence laws consistent with known marginals plus pairwise information.
- The decisive negative result (marginals only implies a dictator) comes from economics (de Oliveira–Ishii–Lin; Arieli et al.), not from jury-theorem papers.

### Cited Findings
- Karge & Rudolph, "The More the Worst-Case-Merrier: A generalized Condorcet jury theorem for belief fusion," KR 2022 proceedings (19th KR, pp. 205–214). [FT this session] — [PDF](https://iccl.inf.tu-dresden.de/w/images/4/4c/Kr2022-0021-karge-rudolph.pdf)
  - The worst case is over a family of joint laws satisfying the paper's assumptions, but agents' approvals are **conditionally independent given the state**. Dependence is allowed only between the scores of different choices.
  - Asymptotic convergence to 1, a finite-n bound, and failure of the non-asymptotic CJT.
- Karge et al., "To lead or to be led" (AAMAS 2024): opinion-leader dependence. Karge, "Bridging dependence: from opinion leader influence to standard correlation measures," arXiv 2609.00010 (29 July 2026) translates the opinion-leader model into covariance and correlation; the abstract mentions no worst-case bounds or weights. [ABS] — [AAMAS24 PDF](https://iccl.inf.tu-dresden.de/w/images/e/e3/AAMAS24_OL.pdf); [arXiv](https://arxiv.org/abs/2609.00010)
- Althöfer & Thiele, *Theory and Decision* 81(1) (2016): tight bounds with correlated couples (see Q1). [ABS] — [IDEAS](https://ideas.repec.org/a/kap/theord/v81y2016i1d10.1007_s11238-015-9521-0.html)
- Baharad, Nitzan & Segal-Halevi, "On the degree of generalizability of Condorcet jury theorem," *Theory and Decision* 99(4):781–800 (2025), DOI 10.1007/s11238-025-10076-6. Monte Carlo study of heterogeneous competence; no dependence and no worst case. [ABS] — [Springer](https://link.springer.com/article/10.1007/s11238-025-10076-6)
- Nakada, Nitzan & Ui, "Robust voting under uncertainty," arXiv 2507.22655 (2025). Uncertainty is over *preference* profiles; when the uncertainty set contains all distributions, the robust rules are exactly tie-free weighted majority rules. The criterion is Pareto-undominance, not competence. [PRIOR-NOTE] — [arXiv](https://arxiv.org/html/2507.22655)
- Ambiguity-averse strategic voting: Ellis, "Condorcet meets Ellsberg," *Theoretical Economics* (2016), and the AUT working paper "Feddersen and Pesendorfer meet Ellsberg" (2019; later versions to 2024). Strategic voters with maxmin expected utility over ambiguous signal likelihoods; equilibrium voting, not weight design. [Search listing only, not read] — [LSE eprint](https://eprints.lse.ac.uk/64598/); [AUT WP](https://www.aut.ac.nz/__data/assets/pdf_file/0011/308558/1907-Feddersen-and-Pesendorfer-meet-Ellsberg.pdf)
- Kaniovski & Zaigraev, *Theory and Decision* 71(4):439–459 (2011). A "correlation-robust" quota for homogeneous exchangeable juries; this is invariance, not max–min. [PRIOR-NOTE] — [IDEAS](https://ideas.repec.org/a/kap/theord/v71y2011i4p439-459.html)
- Bai, Wang, Muller, Cheng & Chen, "Stability of weighted majority voting under estimated weights," AAMAS 2024 (arXiv 2207.06118). Log-odds weights with estimated competences, **assuming independence**; dependent sources are left as future work. [FT this session] — [arXiv](https://arxiv.org/pdf/2207.06118)

### Inferences
- The social-choice literature gives (i) the K=3 / Kaniovski benchmark numbers for validation and (ii) the "dependence-invariant rule" notion (Kaniovski–Zaigraev quota; the notes' Δ(h)=0 rules). It does not threaten the novelty of max–min weights under pairwise information.
- Bai et al. (AAMAS 2024) is a useful foil: robustness to *estimation* error in competences under independence. The project's setting is robustness to *dependence*, and a paper could combine both via confidence sets on p̂.

### Gaps
- Not read: Ellis (2016) and the AUT working paper; Berg (1994); Shapley & Grofman (1984); Nitzan & Paroush (1984, "total interdependence"); Ladha (1995); Pivato (2017). Berg 1994 is the most likely to contain overlapping "weighted rules under dependent voting" results.

---

## Q5. Probability/OR side: Fréchet bounds for weighted Bernoulli sums and threshold events, DRO chance constraints with marginal ambiguity, and complexity

### Takeaway
- Building blocks for (a) and (b) are classical:
  - Hailperin LP;
  - Rüger/Morgenstern closed form for k-of-n given marginals, with a constructive proof by Salako 2025;
  - Prékopa / Boros–Prékopa binomial-moment closed forms;
  - tree-structured pairwise LPs (Natarajan–Padmanabhan);
  - closed forms under pairwise independence (Ramachandra–Natarajan);
  - pseudo-polynomial compact LPs for weighted sums (Padmanabhan et al.);
  - NP-hardness landmarks.
- **Not found anywhere:**
  - sharp bounds for a *weighted* threshold (weighted-majority) event given individual marginals **plus pairwise** probabilities;
  - any DRO formulation that *optimises the decision rule (weights)* against a Fréchet + pairwise ambiguity set for a threshold event.
- DRO chance-constraint work optimises continuous decisions under moment or Wasserstein sets, not Boolean threshold rules.

### Cited Findings
- Natarajan & Padmanabhan, "Tree bounds for sums of Bernoulli random variables: a linear optimization approach," *INFORMS J. Optimization* (DOI 10.1287/ijoo.2019.0038; arXiv 1910.06321).
  - Tight bounds on P(ΣZ_i ≥ k) given p_i and p_ij on a tree, as an O(n²)-variable LP.
  - Covers weighted sums of count probabilities, not Σc_iZ_i with variable weights.
  - [PRIOR-NOTE] — [arXiv](https://arxiv.org/html/1910.06321)
- Ramachandra & Natarajan, "Tight probability bounds with pairwise independence," arXiv 2006.00516 (v6, 23 Nov 2022). [FT abstract page this session] — [arXiv](https://arxiv.org/abs/2006.00516v6)
  - Closed-form tight union bound under pairwise independence; the Boole/tight ratio is ≤ 4/3.
  - For k ≥ 2, new non-tight upper bounds, tight for identical marginals.
  - The arXiv page lists no journal reference; an SUTD repository lists a journal article version that I did not resolve.
- Salako, "Constructive proofs of generalized Boole–Fréchet bounds: a DP approach," arXiv 2512.09161 (Dec 2025). Weighted and threshold events are not treated. [PRIOR-NOTE] — [arXiv](https://arxiv.org/html/2512.09161v1)
- Salako & Muhammad, "Fixed-point characterisations of extremal distributions under partial distributional constraints," arXiv 2608.04315 (5 Aug 2026). A general robust-inference framework for ratio-of-expectation objectives; not specific to Bernoulli thresholds. [ABS] — [arXiv](https://arxiv.org/abs/2608.04315)
- Padmanabhan, Ahipasaoglu, Ramachandra & Natarajan, arXiv 2109.01591. Hardness and a pseudo-polynomial algorithm, as in Q1. [PRIOR-NOTE] — [arXiv](https://arxiv.org/pdf/2109.01591)
- Fu & Li, "Fourth-order marginal moment model: reformulations and applications" (Operations Research, forthcoming per the template; year not printed).
  - Marginal moment models give SDP/SOCP reformulations for worst-case expected objectives of 0-1 problems.
  - Cites NP-hardness of persistency with two-valued coefficients (Bertsimas et al. 2006).
  - No Bernoulli threshold or voting applications.
  - [FT this session] — [optimization-online](https://optimization-online.org/wp-content/uploads/2024/02/4thMomentModel_-7.pdf)
- Embrechts–Puccetti (2006), Puccetti–Rüschendorf (2012/13), Bernard–Rüschendorf–Vanduffel (JRI 2017; variance constraints), Bartl et al. (SIAM, OT with dependence information), Dhara–Das–Natarajan (*IJOC* 2021). Risk-aggregation Fréchet bounds with extra information; continuous or VaR focus. [PRIOR-NOTE]
- Agrawal, Ding, Saberi & Ye (Operations Research), the correlation gap; Dughmi (survey): convex and concave closures. Evaluating the concave closure is NP-hard. [PRIOR-NOTE] — [ar5iv](https://ar5iv.arxiv.org/html/0902.1792)
- Huber & Marić (arXiv 1706.06182); Chaganty & Joe (*Biometrika* 2006); Embrechts, Hofert & Wang (*AAP*). Bernoulli correlation compatibility is a cut-polytope question; CUT_6 has 368 facets, which is relevant to K=5 with pairwise data. [PRIOR-NOTE]
- DRO joint chance constraints (e.g. Xie & Ahmed; Hanasusanto et al.; Bonferroni approximations). These use moment or Wasserstein sets with continuous decisions; I saw none with Fréchet (fixed-marginal) sets over Boolean threshold rules. [Search listing only] — [optimization-online](https://optimization-online.org/wp-content/uploads/2017/02/5860.pdf); [SIAM](https://epubs.siam.org/doi/10.1137/16M1094725)

### Inferences
- The "covering formula" in the notes (worst-case accuracy of a monotone rule given marginals = 1 − min(1, τ*_q), with τ*_q a fractional covering number) is almost certainly folklore in the Hailperin/Boolean-probability literature. It is a monotone-event LP whose dual is a fractional covering. Present it as a lemma with citations to Hailperin, Rüger and Hart–Weiss, not as a contribution.
- The OR side is where a theory reviewer would expect a K-general result. Two candidates that appear open:
  - a pseudo-polynomial DP for the worst-case accuracy of a weighted vote given a **tree** of pairwise probabilities, extending Natarajan–Padmanabhan from counts to weighted thresholds;
  - hardness of the outer max over weights (e.g. Σ₂-type or NP-hard even for fixed small K with integer weights).

### Gaps
- I could not read Fu & Li's publication year, the full text of Boros–Prékopa 1989, or the Embrechts–Puccetti dual formulas.
- A search on "Extended convolution bounds on the Fréchet problem" (arXiv 2511.21929, risk aggregation) returned an unreadable PDF.

---

## Q6. Overall verdict: per-candidate novelty risk and what would make it NeurIPS-worthy

### Takeaway
| Candidate | Closest prior work | What is new | Novelty risk |
|---|---|---|---|
| a1: marginals only → single best classifier | de Oliveira–Ishii–Lin Thm 1; Mazzetto Prop 3; Arieli et al. 2025; B&F 2015 | Nothing (only a one-line covering proof) | **HIGH** (known) |
| a2: K=3 max–min rule with pairwise data | Mazzetto AISTATS 2021 Prop 5 + subset selection | T-form, symmetric phase boundary, signed rules | **HIGH** |
| a3: K=5 phase diagram / structure | Mazzetto PGMV (unweighted odd subsets, heuristics); de Oliveira §7.1.2 open | Weighted and signed classes; a *proved* "committee size ≤ 3" theorem or counterexample | **MEDIUM–HIGH** |
| a4: general-K complexity | Kreinovich–Ferson; Padmanabhan et al.; Pitowsky; Kaski et al. ESA 2026 | Outer max–min hardness and approximation; tree-structured DP for weighted thresholds | **MEDIUM** |
| b1: even/odd identification theorem | Platanios 2014; spectral sign ambiguity; B&F; Corrada-Emmanuel 2025 | Exact identified set E(m,p); only single-classifier accuracies identified; ≤K fractional pairs | **MEDIUM** (new statement, known ingredients) |
| b2: sharp intervals for weighted votes | Mazzetto 2021 LP (unweighted); Maia Polo NeurIPS 2024 (classifier metrics); Feofanov JMLR 2024 | Weighted/signed votes, pair-splitting LP plus dual certificate, finite-sample sets | **MEDIUM** |
| b3: sharp vs PAC-Bayes bounds | **Liu 2026 Thm 4.3 (same formula)**; Boros–Prékopa 1989; Germain 2015 Prop 13; Theisen 2023 lower bound | L(MV) ≤ E[L]+E[D] in PAC-Bayes form; proof that C-bound, tandem and first-order bounds are never sharp for finite K; weighted V_ρ; per-classifier (p,r) bounds | **HIGH** as a new inequality, **MEDIUM** as a sharpness framework |
| b4: robust or minimax-regret weights | B&F 2015 (randomised, transductive); Mazzetto | Deterministic weighted votes, regret criterion, inductive | **MEDIUM** |
| c: LLM juries / voting application | Kohli 2026; J. Chen 2026; Ai ICML 2026; CARE ICML 2026; Liu 2026 | Assumption-free certified intervals and max–min weights for heterogeneous panels | **LOW–MEDIUM** (as application of a+b); **HIGH** if only "LLMs are correlated" or "voting is non-monotone" |
| c: weak supervision under misspecification | Maia Polo 2024; Mazzetto 2021; Cachay 2021 | Weighted votes with abstention under Fréchet/pairwise sets | **MEDIUM** |

### Cited Findings
- Each row's evidence is cited in Q1–Q5 above. Key decisive citations:
  - [de Oliveira–Ishii–Lin](https://arxiv.org/html/2106.00088)
  - [Mazzetto et al. 2021](https://par.nsf.gov//servlets/purl/10273522)
  - [Liu 2605.03379](https://arxiv.org/pdf/2605.03379)
  - [Theisen et al. 2023 supplement](https://proceedings.nips.cc/paper_files/paper/2023/file/30b6fa308e62ed52180c31ae3ba6bb0a-Supplemental-Conference.pdf)
  - [Lacasse et al. 2006](http://www.pascalgermain.info/publis/nips06a.pdf)
  - [Germain et al. 2015](https://ar5iv.arxiv.org/html/1503.08329)
  - [Maia Polo et al. 2024](https://arxiv.org/pdf/2312.04601)
  - [Corrada-Emmanuel 2510.00821](https://arxiv.org/html/2510.00821v1)
  - [Kohli 2026](https://arxiv.org/html/2605.29800v1)
  - [J. Chen 2026](https://arxiv.org/abs/2606.27288)

### Inferences
**Recommended positioning: merge b1 + b2 + b4, with a3 as the structural result and c as the testbed.** A working title: "Certified ensembles under unknown dependence: identification, sharp accuracy intervals and max–min weights from few labels and many unlabeled agreements."

**Theory a NeurIPS reviewer would need, ranked by value against novelty risk:**
1. **Identification theorem (b1).** E(m,p) is the identified set; only single-classifier accuracies are identified; with Z ⟂ Y and η ≠ ½, π is point-identified. State the relation to Platanios and B&F explicitly.
2. **Structure theorem for the max–min rule over the pairwise set.** Prove (or refute with an explicit family) that the max–min rule is attained by a committee of size ≤ 3, or characterise when larger or weighted committees win, e.g. under negative dependence or the "anti" family.
   - Proved: a strong result that also explains why Mazzetto's size-3 heuristic works.
   - Refuted: the counterexamples (weighted (2,1,1,1,0)-type winners) are the novelty over Mazzetto.
3. **Pair-splitting LP + dual certificate (b2)** with finite-sample confidence sets (Clopper–Pearson or PAC-Bayes-kl on p̂; concentration on agreement), and a label-budget analysis. The notes' empirical finding that labeled uncertainty dominates is a strong practical message.
4. **Sharpness section (b3),** framed honestly:
   - cite Liu 2026 Thm 4.3 and Boros–Prékopa for the moment-problem solutions;
   - the contributions are (i) the corollary L(MV_ρ) ≤ E_ρ[L] + E_{ρ²}[D], the mirror of Theisen's lower bound; (ii) a theorem that the C-bound, tandem and first-order bounds are strictly non-sharp for every finite ensemble; (iii) per-classifier bounds that dominate averaged-moment bounds.
5. **Complexity (a4).** NP-hardness of the outer max–min, plus a tractable tree-structured DP. Optional.

**Experiments needed:**
- (i) Synthetic sharpness gaps versus FO, TND, C-bound, CCPBB, Theisen and Kim bounds.
- (ii) Classical ensembles where errors are *complementary*: multi-view or feature-subset ensembles, plus covariate shift with unlabeled target data. The notes' scikit-learn study showed that same-feature ensembles collapse to a single classifier.
- (iii) LLM judge panels on binary or pairwise-preference tasks (RewardBench, UltraFeedback pairs, ChaosNLI-binarised), deliberately including same-family redundant judges. Baselines: majority vote, log-odds weights, Dawid–Skene, Ai et al. OW/ISP, Ising, CARE and Mazzetto PGMV. Report both certified worst-case and realised accuracy.
- (iv) WRENCH, with an abstention extension.

**Must-cite and must-differentiate list:**
- de Oliveira–Ishii–Lin (incl. §7.1.2);
- Mazzetto et al. AISTATS 2021 (Prop 3, Prop 5, PGMV);
- Balsubramani–Freund COLT 2015 / NeurIPS 2016;
- Platanios UAI 2014;
- Maia Polo NeurIPS 2024;
- Liu arXiv 2605.03379;
- Germain JMLR 2015; Lacasse NIPS 2006; Masegosa NeurIPS 2020; Wu NeurIPS 2021;
- Theisen NeurIPS 2023; Kim NeurIPS 2024;
- Feofanov JMLR 2024;
- Corrada-Emmanuel 2025;
- Kohli 2026; J. Chen 2026;
- Arieli et al. AEJ:Micro 2025;
- Boros–Prékopa 1989; Natarajan–Padmanabhan.

### Gaps
- No direct OpenReview full-text search for NeurIPS 2025/2026 or ICLR 2026 submissions was possible, so recent accepted-but-not-on-arXiv overlaps cannot be ruled out.
- Older items flagged as possibly overlapping and not read: Narasimhamurthy TPAMI 2005; Berg 1994; Lam & Suen 1997; Boros–Prékopa 1989 full text.
