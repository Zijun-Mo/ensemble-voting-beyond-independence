# Novelty check: recent venues (2024-2026) for worst-case aggregation of dependent binary classifiers / LLM judges

Date of search: 2026-10-13. Tools: WebSearch and WebFetch only. Read-depth tag: all items are [ABS] (abstract or partial HTML excerpt via summarizer; no full PDF read, so no theorem statements were verified directly).

## Ranked list (at most 15, by overlap with our question)

1. **Dependence-Aware Label Aggregation for LLM-as-a-Judge via Ising Models**
   Krishnakumar Balasubramanian, Aleksandr Podkopaev, Shiva Prasad Kasiviswanathan. ICML 2026 (PMLR v306). arXiv 2601.22336 (Jan 2026).
   Summary: Models K binary judges with class-conditional Ising (pairwise MRF) dependence and latent-factor dependence. For class-independent couplings, the Bayes log-odds reduce to a linear weighted vote with correlation-adjusted weights. Shows separation results: conditional-independence aggregators stay suboptimal as K grows under latent factors.
   Overlap: **HIGH**. Same object (pairwise dependence, weighted vote, LLM judges). Difference: Bayes-optimal under a known parametric Ising model, not worst-case over unknown joint distributions; no odd-cycle or co-error-graph structure seen; no Cantelli/C-bound.
   [ABS]

2. **CARE: Confounder-Aware Aggregation for Reliable LLM Evaluation**
   Jitian Zhao, Changho Shin, Tzu-Heng Huang, Satya Sai Srinath Namburi GNVV, Frederic Sala. ICML 2026 (PMLR v306). arXiv 2603.00039 (v2 Sep 2026). Workshop precursor: "From Many Voices to One: A Statistically Principled Aggregation of LLM Judges", NeurIPS 2025 workshop "Reliable ML from Unreliable Data" (listed as arXiv 2512.03019 is a different paper; the workshop version was not separately fetched).
   Summary: Models judge scores as a latent true-quality signal plus k shared confounders; sparse-plus-low-rank precision estimation (CARE-SVD, Gaussian) and tensor CP decomposition (CARE-Tensor, binary latents). Aggregates with identified loadings. Claims identifiability and finite-sample recovery; reports up to 26.8% error reduction on 12 benchmarks.
   Overlap: **MEDIUM-HIGH** on the "dependence-aware judge aggregation" framing. Not worst-case; parametric (Gaussian / Gaussian mixture); no pairwise-disagreement bounds; no weighted-majority optimality claim.
   [ABS]

3. **Nine Judges, Two Effective Votes: Correlated Errors Undermine LLM Evaluation Panels**
   Guneet Kohli (single author; Apple research blog also lists it). arXiv 2605.29800 (May 2026). No venue stated.
   Summary: 9-model LLM-judge panel on NLI (ChaosNLI) and RewardBench. Mean pairwise phi of errors; Kish effective sample size n_eff = k/(1+(k-1)phi); panel worth about 2 independent votes; Condorcet null via Monte Carlo confusion matrices. States that Dawid-Skene and accuracy-weighted voting close at most 11% of the gap. Does not derive a worst-case or optimal-rule theorem; argues empirically that no weighting extracts a third independent vote.
   Overlap: **MEDIUM** on motivation and metric (pairwise error correlation → lost votes). No robust or structural aggregation result. Useful as motivation citation.
   [ABS]

4. **Robust Decision Aggregation with Adversarial Experts**
   Yongkang Guo, Yuqing Kong. arXiv 2403.08222 (Mar 2024; v2 Feb 2025). No venue stated in the page.
   Summary: Minimax-regret aggregation of binary expert reports when some experts are adversarial, with marginally symmetric experts and unknown information structure. Truncated mean is optimal under L1; for hard decisions, a randomized truncated mean is optimal for both L1 and L2. Regret depends on the adversary ratio.
   Overlap: **MEDIUM-HIGH** on the method class (worst-case/minimax over unknown structure, randomized aggregator optimal for decisions). Differences: adversaries rather than pairwise dependence; unweighted truncated mean, not weighted majority; no odd-cycle structure.
   [ABS]

5. **No-Knowledge Alarms for Misaligned LLMs-as-Judges**
   Andrés Corrada-Emmanuel. arXiv 2509.08593 (Sep 2025). Note: this is an earlier arXiv entry than the already-known 2510.00821 ("Logical Consistency Between Disagreeing Experts and Its Role in AI Safety"); the two appear to be the same line of work. Please verify whether 2509.08593 is already in your list.
   Summary: Uses logical consistency between disagreeing judges (e.g., if A and B disagree about C's correctness, not both can be fully correct). Encoded as an LP over integer response counts. Produces "no-knowledge alarms" with no false positives when a judge violates a user-specified accuracy requirement.
   Overlap: **MEDIUM-HIGH** on the bounds-from-disagreement-counts mechanism. Detection rather than aggregation; no weighted vote; no worst-case aggregation rule.
   [ABS]

6. **Partial Identification from LLM Prompts**
   Xiaohong Chen, Ashesh Rambachan, Elie Tamer. arXiv 2606.15031 (June 2026; econ.EM). No venue.
   Summary: Prevalence of a latent binary label when LLM classifiers have arbitrarily dependent errors. Shows the prevalence is unidentified without separating restrictions; identification comes from externally calibrated scores/events. Sharp bounds; key point that which models agree across prompts carries identifying information that a vote count discards.
   Overlap: **MEDIUM**. Worst-case / sharp bounds under arbitrary dependence on LLM classifier outputs, and agreement patterns are informative. It targets prevalence, not ensemble accuracy or the choice of aggregation rule.
   [ABS]

7. **Unsupervised Estimation of Ensemble Accuracy**
   Simi Haber, Yonatan Wexler. arXiv 2311.10940 (Nov/Dec 2023; NeurIPS 2023 workshop / "Practical estimation of ensemble accuracy", NeurIPS 2023). Outside the 2024-2026 window; included as a borderline item.
   Summary: Label-free combinatorial bound on the number of ensemble mistakes, approximated in linear time in samples. Demonstrated on face recognition with noisy labels.
   Overlap: **LOW-MEDIUM** (label-free ensemble accuracy estimation). No pairwise bound or robust aggregation.
   [ABS]

8. **RoPoLL: Robust Panel of LLM Judges**
   Anish Acharya, Kris W. Pan, Brian Verkhovsky (AWS). arXiv 2606.30931 (June 2026). Accepted to ICML 2026 Workshop on Statistical Frameworks for Uncertainty in Agentic Systems.
   Summary: Huber contamination model for LLM jury (PoLL). Shows arithmetic-mean jury bias is unbounded under contamination; proposes geometric median; finite-sample bound and minimax lower bound (rate sigma*sqrt(d/N)). Continuous scores, not binary classifiers. Lemma 3 extends to equicorrelated juries.
   Overlap: **LOW-MEDIUM**. Worst-case/minimax aggregation of judges, but over contamination of continuous scores; not dependence bounds, not weights, not graph structure.
   [ABS]

9. **Don't Always Pick the Highest-Performing Model: An Information Theoretic View of LLM Ensemble Selection**
   Yigit Turkmen, Baturalp Buyukates, Melih Bastopcu. arXiv 2602.08003 (Feb 2026).
   Summary: Models correlated LLM errors with a Gaussian copula; derives an information-theoretic error floor; greedy mutual-information selection under a query budget. Average-case, not worst-case.
   Overlap: **LOW-MEDIUM** (correlated errors of ensemble members; selection, not aggregation rule).
   [ABS]

10. **The Majority Vote Paradigm Shift: When Popular Meets Optimal**
    Antonio Purificato, Maria Sofia Bucarelli, Anil Kumar Nelakanti, Andrea Bacciu, Fabrizio Silvestri, Amin Mantrach. AISTATS 2026 (PMLR v300). arXiv 2502.12581 (v4 Feb 2026).
    Summary: Conditions under which unweighted majority vote achieves the optimal lower bound on label estimation error; tolerable noise levels as a function of class distribution.
    Overlap: **LOW-MEDIUM**. Optimality of MV for aggregation; no dependence, no weights, no worst-case over joint distributions.
    [ABS]

11. **Are Diversity Metrics Measuring Diversity? A Capability-Controlled Audit of Majority-Vote Gain in LLM Ensembles**
    Donghwan Kim. arXiv 2607.20768 (July 2026).
    Summary: Empirical audit over 31,900 LLM subsets (MMLU-Pro). Oracle gain is positive, but simple voting beats the best member rarely. After capability control, the only stable signal is a modest pairwise co-failure association (more shared errors, lower gain).
    Overlap: **LOW-MEDIUM** (empirical support for pairwise co-error as the relevant statistic). No theory.
    [ABS]

12. **Quantifying Diversity of Thought: A Predictive Law of Weighted LLM Ensemble Lift**
    Junade Ali. arXiv 2607.17384 (July 2026; v3 Sept 2026).
    Summary: Decomposes ensemble lift into rescue and damage masses; heuristic based on accuracy-adjusted pairwise correctness correlation (phi_adj), accuracy gap, and pair collective accuracy. Empirical calibration on SuperGPQA, tested on GPQA Diamond and a forensic benchmark.
    Overlap: **LOW-MEDIUM** (pairwise correctness correlation predicts weighted-vote lift). Heuristic, not worst-case; no structural characterization.
    [ABS]

13. **Generalizing Condorcet's Jury Theorem to Social Networks**
    Dan Braha, Marcus A. M. de Aguiar. arXiv 2510.16808 (Oct 2025). Physics/social network venue.
    Summary: Jury theorem under network imitation dependence with zealots; exact accuracy limits for majority and plurality.
    Overlap: **LOW**. Dependence from a graph, but unweighted rules and a specific dynamic model; no worst-case over structures.
    [ABS]

14. **LLMs Judging LLMs: A Simplex Perspective** (arXiv 2505.21972; the link we were asked to fetch)
    Patrick Vossler, Fan Xia, Yifan Mai, Adarsh Subbaswamy, Jean Feng. arXiv 2025 (v3 Apr 2026). No venue.
    Summary: Geometric/Bayesian priors on judge quality for judge-only ranking of LLM outputs; identifiability conditions for rankings.
    Overlap: **LOW**. Not about aggregation of judges or dependence. Note: 2505.21972 is this paper, not a paper on ensemble accuracy bounds.
    [ABS]

15. **Epistemic Filtering and Collective Hallucination: A Jury Theorem for Confidence-Calibrated Agents**
    Jonas Karge. arXiv 2602.22413 (Feb 2026; v2 Apr 2026).
    Summary: Abstention via confidence gating; non-asymptotic lower bound on group success under fixed competences. Voters independent (dependence not modeled).
    Overlap: **LOW**. Jury-theorem bound, but no dependence or worst-case.
    [ABS]

### Checked and judged no or negligible overlap (not in top 15)
- Majority-of-Three is Optimal (Rawal & Zhivotovskiy, arXiv 2606.13614, June 2026): realizable PAC, independent classifiers; no dependence or weights. Overlap: none.
- Majority-of-Three: The Simplest Optimal Learner? (Aden-Ali et al., COLT 2024, arXiv 2403.08831): same; none.
- Multi-class Probabilistic Bounds for Majority Vote Classifiers with Partially Labeled Data (Feofanov, Devijver, Amini; ICML 2024 / JMLR 2024): majority-vote risk bounds for self-training; no dependence model; none.
- Validating LLM-as-a-Judge Systems under Rating Indeterminacy (Guerdan et al., arXiv 2503.05965): human-label aggregation; none.
- Distribution-Calibrated Inference Time Compute for Thinking LLM-as-a-Judge (Dadkhahi et al., arXiv 2512.03019): repeated samples from one judge; none.
- A Study of LLM Judge Capability-Dependent Biases and Multi-Judge Ensemble for Bias Calibration (NeurIPS 2026 listing): not fetched; title suggests bias calibration, not dependence; low.
- Bridging Dependence: From Opinion Leader Influence to Standard Correlation Measures (arXiv 2609.00010): surfaced in search but NOT fetched; unverified, title suggests a social-influence model; low (unverified).

## Searches that failed or were inconclusive
- **OpenReview:** not directly searchable through WebSearch. Queries for NeurIPS 2025 / ICLR 2026 / NeurIPS 2026 returned almost nothing from OpenReview itself. OpenReview coverage is therefore a real gap; NeurIPS 2026 submissions are probably not public.
- **COLT 2025 / 2026:** no relevant hits (results were 2015-2024 items, e.g., Balsubramani & Freund COLT 2015, already known). Result 2601.02193 ("Learning with Monotone Adversarial Corruptions") surfaced but was not fetched.
- **AISTATS 2025 / UAI 2025 and proceedings.mlr.press browsing:** not systematically browsed. Only AISTATS 2026 (Purificato et al.) and ICML 2026 (Ising, CARE) items were hit by search.
- **Odd cycles:** no recent paper found combining ensemble voting with odd cycles, co-error graphs, or bipartiteness. Search hits were all social choice or old voting papers. Strong null signal, but not exhaustive.
- **C-bound / Cantelli, recent:** only 2014-2015 C-bound generalizations (Lacasse et al. line) and a 2021 NeurIPS paper. No 2024-2026 C-bound paper with dependence-based worst-case analysis was found. No recent "half Cantelli" hit.
- **Fréchet bounds for ensembles:** only 2016-2023 items (e.g., "Estimating Fréchet bounds for validating" arXiv 2312.04601, 2023, not fetched in depth). Nothing 2024-2026 applying Fréchet bounds to classifier co-error.
- **"Triplet" / label-free accuracy with dependence (Jaffe et al. line):** no 2024-2026 extension found beyond the known 2014-2015 works.
- **Known list items** (Mazzetto, de Oliveira-Ishii-Lin, Arieli et al., Liu 2605.03379, Theisen et al., Kim et al. 2024 "How many classifiers", Maia Polo et al. NeurIPS 2024, Corrada-Emmanuel 2510.00821, Kim et al. ICML 2025 "Correlated Errors in LLMs") were not re-reported.

## Verdict: has anyone in 2024-2026 posed or solved our question?

**No exact match found.** Our specific combination is not in the recent literature as far as these searches show:
- unknown joint error distribution (distribution-free / adversarial),
- known marginal accuracies plus pairwise co-error upper bounds estimated from unlabeled disagreement,
- a worst-case-optimal aggregation rule (weighted majority, subset majority, randomized committee),
- structural characterization via odd cycles in a low-co-error graph (bipartite gives no gain), and a randomized half-Cantelli / half-C-bound rule for large K.

**Closest threats (in order of concern):**
1. **Balasubramanian et al., ICML 2026 (arXiv 2601.22336).** Weighted vote with correlation-adjusted weights under pairwise Ising dependence, plus separation results. This is the main threat to any "weighted vote under dependence" framing. Our differentiator must be the worst-case / distribution-free setting and the odd-cycle structure; the Ising paper is Bayes-optimal under a parametric model.
2. **Guo & Kong (arXiv 2403.08222).** Randomized truncated mean is minimax-regret optimal for hard decisions with adversaries. Closest precedent for a randomized worst-case aggregator; the difference is no pairwise dependence and no weights.
3. **Chen, Rambachan & Tamer (arXiv 2606.15031).** Sharp worst-case bounds under arbitrary dependence for LLM classifiers. Different target (prevalence) but the same worst-case-under-arbitrary-dependence spirit; cite as a related framework.
4. **CARE (ICML 2026) and Nine Judges (arXiv 2605.29800).** Motivation and modeling overlap; neither does worst-case rules or pairwise-bound optimization.
5. **Corrada-Emmanuel, arXiv 2509.08593.** Disagreement-count LP bounds; likely already known via 2510.00821, but the 2509 version may not be.

**Bottom line:** The field has converged on "dependence matters" for LLM judges (CARE, Nine Judges, Ising, Kim ICML 2025), and on worst-case / minimax thinking in adjacent settings (Guo & Kong, RoPoLL, Chen et al.). The **worst-case aggregation rule over pairwise co-error bounds with odd-cycle characterization and a randomized half-Cantelli / half-C-bound rule** does not appear in 2024-2026 venues from these searches. Coverage gaps (OpenReview, COLT/AISTATS/UAI 2025, NeurIPS 2026) mean this is a "not found" result, not a proof of novelty. Recommend a targeted check of OpenReview ICLR 2026 and NeurIPS 2025 submissions by title before claiming novelty.
