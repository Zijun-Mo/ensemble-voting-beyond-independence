# Novelty check 2: minimax aggregation under pairwise co-error bounds

Scope: check whether Claims A1, A2, A4 (self-dual minimax aggregation with known q_i and pairwise co-error upper bounds b_ij) are already in the literature.

Source marks:
- [FT] full PDF was read through the WebFetch extraction (a model summary of the whole PDF, not line-by-line verification of every theorem).
- [FT-partial] only part of the PDF or page was reached.
- [ABS] abstract or landing page only.
- [NA] not accessible or not read.

---

## 1. Balsubramani and Freund line (first-order constraints)

**1a. Balsubramani & Freund, "Optimally Combining Classifiers Using Unlabeled Data", COLT 2015. arXiv:1503.01811. PMLR v40.** [FT]
- Setting: p binary classifiers, n unlabeled test points, prediction matrix F. Only label information is per-classifier correlation lower bounds b_i (from labeled data). Constraint: (1/n) F z >= b, where z in [-1,1]^n is the adversary's labeling. No pairwise or higher-order constraints.
- Main result (Thm 2 in the extraction): minimax value V = -gamma(sigma*), with slack gamma(sigma) = (1/n) sum_j [|x_j' sigma| - 1]_+ - b' sigma, minimized over sigma >= 0. Gamma is convex, so sigma* is an LP or SGD solution.
- Aggregator: nonnegative weighted vote x' sigma*, clipped to [-1,1]. Predictions are randomized (hedged) where |margin| < 1. With F unknown, the solution reduces to the best single classifier, max_i b_i.
- Implies A1? Partially. Their "linear weighted vote with nonnegative weights" is the first-order special case (gamma = 0) of A1. However their rule is clipped and randomized, not deterministic self-dual. The sign reduction for deterministic self-dual rules is not stated.
- Implies A2? Partially. Their nonnegativity of sigma comes from correlation constraints b >= 0. Our condition sum q_i <= 1 is not used or discussed.
- Implies A4? Partially. The "best single classifier" branch matches their max_i b_i ERM reduction. Their paper has no graph bound and no q_min/B_graph dichotomy.

**1b. Balsubramani & Freund, "Scalable Semi-Supervised Aggregation of Classifiers", NeurIPS 2015. arXiv:1506.05790.** [FT]
- Same game as 1a with a scalable SGD / minibatch method. The extraction says it "builds directly on" COLT 2015 and does not contrast with it. Constraints are per-classifier only. Output is a clipped weighted vote.
- NeurIPS proceedings abstract page (papers.nips.cc ce78d1da...) [ABS]: no constraint or form details beyond "no assumptions on structure or origin."
- Implies A1, A2, A4: same as 1a, partially for A1 and A2 only.

**1c. Balsubramani & Freund, "Optimal Binary Classifier Aggregation for General Losses", NeurIPS 2016. arXiv:1510.00452.** [FT]
- Same first-order setting with general partial losses l+ and l-. Theorem 4 gives V = (1/2) min over sigma >= 0 of [-b'sigma + (1/n) sum_j Psi(x_j' sigma)], with a loss-dependent potential Psi. Predictions are a clipped sigmoid-like function of the weighted margin, and they are randomized in [-1,1]. The 0-1 loss case gives clipping. The extraction says there are no pairwise terms. Nonnegativity holds for sigma in the main theorem. Prop. 7 and Thm. 11 drop the nonnegativity with an L1 penalty for two-sided constraints.
- Implies A1? Partially, in the first-order case. Pairwise constraints are not present.
- Implies A2? Partially, with the same caveat as 1a.
- Implies A4? No.

**1d. Balsubramani & Freund, "Muffled Semi-Supervised Learning", arXiv:1605.08833 (2016).** [FT]
- Transductive, first-order, per-classifier constraints. The authors state that second-order methods are an open problem. Output is a clipped nonnegative weighted combination (Theorem 1).
- Implies A1, A2, A4: partially (first-order special case only). Also confirms that pairwise / second-order constraints were left open by this group as of 2016.

**1e. Balsubramani, "Learning to Abstain from Binary Prediction", arXiv:1602.08151 (NeurIPS 2016).** [FT]
- Per-classifier constraints only, with an abstention budget alpha. Optimal aggregator is a weighted vote sgn(x_j' sigma*), abstaining with probability 1 - min(1, |x_j' sigma*|). The extraction says pairwise / joint-distribution constraints are not mentioned.
- Implies A1, A2, A4: no for pairwise terms; partially for the weighted-vote form.

**Summary of the Balsubramani line.** All four papers use first-order constraints (1/n) F z >= b. The minimax rule is a nonnegative clipped weighted vote, and the best-single-classifier fallback appears (max_i b_i). None of them has a pairwise term, a deterministic self-dual restriction, a sum q_i <= 1 condition, or a graph-based value. This is the closest prior work and the most important comparison for A1 and A2. A key difference: their adversary chooses labels z on a fixed empirical pattern set with first-order constraints. Our adversary chooses the joint law of X over {0,1}^K with second-order upper bounds, which is a different LP.

---

## 2. Follow-ups with pairwise / second-order information

**2a. Masegosa, Lorenzen, Igel, Seldin, "Second Order PAC-Bayesian Bounds for the Weighted Majority Vote", arXiv:2007.13532 (NeurIPS 2020).** [FT]
- Uses expected tandem loss over pairs (second-order). Minimizes a high-probability PAC-Bayes bound with softmax (nonnegative, normalized) weights. Not a minimax over joint laws. No worst-case optimization.
- Implies A1, A2, A4: no. Second-order information appears, but as an empirical bound, not a worst-case game.

**2b. Wu, Masegosa, Lorenzen, Igel, Seldin, "Chebyshev-Cantelli PAC-Bayes-Bennett Inequality for the Weighted Majority Vote", arXiv:2106.13624.** [FT]
- Uses first and second moments of the tandem loss. Optimizes weights on the simplex. The offset Chebyshev-Cantelli form is tight over the moment class, but the paper does not frame it as a minimax problem over distributions. Optimal weights are nonnegative normalized.
- Implies A1, A2, A4: no. This is the closest "second-moment + Chebyshev tightness" result, but it is not a worst-case aggregator optimization.

**2c. Germain, Lacasse, Laviolette, Marchand, Roy, "Risk Bounds for the Majority Vote: From a PAC-Bayesian Analysis to a Learning Algorithm" (MinCq), JMLR 16 (2015). arXiv:1503.08329.** [FT-partial] (the extraction did not reach the MinCq section.)
- The C-bound is 1 - mu1^2/mu2 with mu1 = first margin moment and mu2 = second. The extraction states the C-bound is tight given only the first two margin moments (Prop. 13). MinCq minimizes the C-bound over the posterior Q, which "reduces to a simple quadratic program." The weights are a distribution, so they are nonnegative.
- Implies A1: partially. The second-moment structure (pairwise terms) is present. A quadratic-in-weights objective appears, which is similar in spirit to our quadratic certificates. But MinCq optimizes a bound on the risk of a fixed posterior, not a minimax over joint laws, and it has no deterministic self-dual structure.
- Implies A2 or A4: no.

**2d. Jaffe, Nadler, Kluger, "Estimating the Accuracies of Multiple Classifiers Without Labeled Data", AISTATS 2015. arXiv:1407.7644.** [FT]
- Conditional independence (pairwise, triplet-wise, or joint). Estimates sensitivity and specificity from the covariance and 3-way tensor. Weights are log-odds style (i-SML). Gives consistency and rates, not worst-case bounds. Non-identifiability is noted (Theorem 2).
- Implies A1, A2, A4: no.

**2e. Jaffe, Fetaya, Nadler, Jiang, Kluger, "Unsupervised Ensemble Learning with Dependent Classifiers", AISTATS 2016. arXiv:1510.05830.** [FT]
- Dependence is modeled through a latent-cluster structure, not pairwise bounds. The aggregator is nonlinear maximum likelihood. The extraction says there is no minimax or worst-case result, and Lemma 3 is a computational hardness result.
- Implies A1, A2, A4: no. Relevant as a contrast: it handles dependence by modeling structure, while our setting uses worst-case over pairwise upper bounds.

**2f. Zhou, Basu, Mao, Platt, "Learning from the Wisdom of Crowds by Minimax Entropy", NeurIPS 2012.** [FT-partial] (Microsoft Research PDF page extraction; the minimax sets are not fully spelled out.)
- Maximum-entropy over labelings subject to count-based moment constraints (item-level and worker-by-class-pair), solved by Lagrangian dual. This is a minimax over labelings and dual variables, not a worst-case error guarantee. The solution resembles skill-based weighting but is not stated as a weighted vote. No worst-case guarantee.
- Implies A1, A2, A4: no. It is adjacent: second-order count constraints in a minimax-style formulation, but with entropy and not worst-case error.

**2g. Bertsimas & Popescu, "Optimal Inequalities in Probability Theory: A Convex Optimization Approach", MIT Sloan WP 4083 (1999). Later SIAM J. Optim. (2005).** [FT-partial] (PDF extraction of the working paper.)
- General framework: worst-case P(X in S) over distributions with given moments, via LP/SDP duality. The dual minimizes E[g(X)] over polynomial g that dominate the indicator of S. Second-order moments (including E[X_i X_j]) are covered. Does not treat thresholds on sums of indicators or unions explicitly.
- Implies A1: partially, in the sense that the quadratic-certificate duality in A1 is the standard polynomial-certificate duality for moment problems. Strong duality holds in the interior case. A1's specific use (self-dual F, linear reduction via Q(x) - Q(1-x)) is not in this paper.
- Implies A2 or A4: no.

---

## 3. Probability bounds with pairwise (second-order) information on union and threshold events

**3a. Yang, Alajaji, Takahara, "Lower Bounds on the Probability of a Finite Union of Events", arXiv:1401.5543.** [FT]
- LP-based sharp lower bound on P(union) given marginals and pairwise sums. Optimal only within the class of bounds that use P(A_i) and pairwise sums.
- Implies A1: no (union event, not a weighted vote). Relevant as an LP-with-pairwise-moments precedent.

**3b. Ramachandra & Natarajan, "Tight Probability Bounds with Pairwise Independence", arXiv:2006.00516 (v6, April 2022).** [FT]
- Upper and lower bounds on P(at least k of n events) with marginals p and pairwise-independent joint probabilities. The tightest bound is an LP over distributions on {0,1}^n, with a dual that has pairwise and univariate multipliers. For k = 1, a closed form min(sum p_i - p_n sum_{i<n} p_i, 1) matches the Hunter/Worsley spanning-tree bound. The extremal distribution is given. For k >= 2 no general tight bound is claimed. Kwerel (1975) is cited as an aggregated LP that is generally weaker.
- Implies A1: partially. It is a worst-case over joint laws with pairwise terms and an LP dual, which is the same mechanism as A1. But the event is "at least k" (union/intersection), not a self-dual rule. The paper's pairwise setting is exact independence, not upper bounds.
- Implies A4: partially. Its dual and tree structure are the nearest analogue to B_graph. However, the B_graph formula (inf over y of sum b_ij y_ij / g(y), with g the min over S of max of edge sums on the two sides of the cut) is not in this paper. The paper's bounds are over spanning trees for union events, not a min-cut-type quantity for majority events.

**3c. Veneziani, "Graph-Based Upper Bounds for the Probability of the Union of Events", Electronic J. Combinatorics 15 (2008), #R28.** [ABS] (abstract and landing page only.)
- Upper bound on P(union) given individual and pairwise probabilities. The problem is an LP and the bounds are dual basic feasible solutions, with bipartite and threshold-type graphs. The excerpt does not say whether the bounds are optimal over graphs and does not mention cuts.
- Implies A4: partially at most. Graph-based dual characterization for pairwise union bounds. Does not cover self-dual rules or the q_min / B_graph dichotomy.

**3d. Boros, Scozzari, Tardella, Veneziani, "Polynomially computable sharp probability bounds" (Rutgers RUTCOR page, Boros.pdf).** [FT-partial] (web summary of the paper, not the full text.)
- Sharp bounds on P(union) via LPs over joint distributions consistent with pairwise intersection data, with polynomially computable forms. Closed-form S1 and S2 bounds for m = 2 (Chung–Erdős lower bound and Kwerel-type upper bound). Notes that the binomial moment relaxation is not always sharp in terms of pairwise probabilities.
- Implies A4: partially, as a precedent for a sharp, polynomially computable LP-over-joint-laws bound with pairwise data. Not for majority or self-dual rules.

**3e. Kuncheva, "Limits on the majority vote accuracy in classifier fusion", Pattern Analysis and Applications (2003).** [NA] (not read.)
- Note: the URL lucykuncheva.co.uk/papers/cslkif02.pdf that I fetched turned out to be a different Kuncheva paper on diversity measures (Table 3 joint correct/incorrect probabilities). It states that a functional relationship exists between the Q statistic and majority-vote limits, but does not give the bounds. So the 2003 limits paper itself has not been read. Its approach (worst-case or extremal bounds on majority accuracy given pairwise statistics) is the nearest classic precedent for A1/A4 and should be checked in full.

**3f. Kwerel (1975), Hunter (1976), Worsley (1982), Dawson–Sankoff, Boros–Prékopa (1989, closed-form two-sided bounds for at-least-r events via S1, S2).** [NA] (cited via 3b, 3d, and the Boros et al. summary.) These are the classical sharp-bound results for threshold events with binomial moments. Boros–Prékopa is the most relevant for threshold events, but it uses exact or aggregated S2, not per-pair upper bounds b_ij, and its adversary is over pair-sums rather than the full family of pairwise constraints. Needs direct reading to rule out A1 or A4 coverage.

---

## 4. Verdicts

**A1 (threshold sufficiency / linear weighted majority is optimal among self-dual rules, with weights from quadratic certificates): PARTIALLY KNOWN.**
- Closest statements: B-F COLT 2015 Thm 2 and NeurIPS 2016 Thm 4 (first-order case, optimal rule = clipped nonnegative weighted vote, via LP/minimax duality). Bertsimas–Popescu (polynomial-certificate duality for moment problems, second-order moments included). Ramachandra–Natarajan (LP over joint laws with pairwise terms, dual with pairwise multipliers).
- The LP-duality step (worst-case error of a fixed rule = min over quadratic certificates Q >= F) is standard moment-problem duality and is probably not new in itself.
- The specific claim that minimizing over self-dual F reduces to a signed linear vote via Q(x) - Q(1-x) being linear does not appear in any source I read. Balsubramani's clipping/hedging argument is the nearest analogue.
- Caveat to check: B-F's optimal rule is randomized (clipped). A1 claims a deterministic sign rule is optimal. The gap between min over randomized self-dual [0,1]-valued F and min over deterministic self-dual F is not discussed in any source. If there are cases where the fractional relaxation is strictly better, A1 as stated would need the integrality argument.

**A2 (sum q_i <= 1 implies nonnegative weights): NOT FOUND as a stated result.**
- B-F's nonnegativity (sigma >= 0) comes from correlation lower bounds b >= 0 in the first-order setting. Their Prop. 7 / Thm. 11 drop nonnegativity with an L1 penalty for two-sided bounds, which is the closest contrast. No source states a condition in terms of error rates summing to at most 1.
- Partial: the mechanism (nonnegativity of dual multipliers under a sign condition on the constraints) is the same as in B-F.

**A4 (optimal worst-case error = min{q_min, B_graph}; optimal rule is the best single classifier or a weighted-degree vote on a graph): NOT FOUND as a stated theorem. Pieces are known.**
- Known pieces: (i) the best-single-classifier fallback appears in B-F (max_i b_i, which corresponds to q_min in our notation) and is obviously achieved by the dictator rule; (ii) the LP-with-pairwise-terms bounds on union and threshold events have graph-based dual characterizations (Hunter/Worsley spanning tree, Veneziani 2008, Boros et al., Ramachandra–Natarajan 2022).
- Not found: the graph bound B_graph = inf_y sum b_ij y_ij / g(y) with g(y) = min_S max{y(E(S)), y(E(S^c))}, its identification with the self-dual minimax value, and the weighted-degree vote as the optimal rule. The min-over-subsets-of-max-of-cut-sides structure looks like a cut / max-cut-type quantity. I did not find a source that connects a worst-case majority or self-dual error to a cut-type graph quantity.

---

## 5. Items to check before claiming novelty

1. Read Kuncheva (2003), "Limits on the majority vote accuracy in classifier fusion", in full. It is the nearest classic precedent for worst-case majority accuracy given pairwise statistics.
2. Read Boros–Prékopa (1989) and Kwerel (1975) in full. They may already give sharp threshold-event bounds from pairwise sums; if so, A1's LP reduction is not new and only the self-dual / weighted-vote consequence would be.
3. Check randomized vs deterministic: compute, for small K, whether min over randomized self-dual F equals min over deterministic self-dual F under the same constraints. This is the main gap in A1.
4. Check whether B-F's game with n test points, when the empirical pattern distribution is the adversary's variable, reduces to our LP. If it does, that is an easy citation and an honest "special case" statement for A1 and A2.
5. Search Google Scholar "cited by" for Balsubramani–Freund COLT 2015 and for Kuncheva 2003 with the words "pairwise", "worst case", and "self-dual". I did not get live citation lists; the search tool only returned the original papers.
6. Note the limits of this check: Google Scholar citation pages and full texts behind paywalls (Pattern Analysis and Applications, SIAM J. Optim., EJC PDFs) were not directly read. Some sources are verified only through WebFetch extraction.
