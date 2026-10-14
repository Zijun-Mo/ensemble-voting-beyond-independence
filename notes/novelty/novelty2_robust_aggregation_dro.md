# Novelty check 2: robust information aggregation and DRO with pairwise co-error bounds

Date: 2026-10-13
Status key: [FT] full text fetched (the fetch tool returns an extracted/summarized version, so spot-check quotes before citing); [ABS] abstract or metadata only; [NA] not accessible or not fetched.

## 0. Problem being checked

K binary signals (classifier votes) about a binary state. Known: accuracies p_i (errors q_i = 1 - p_i) and upper bounds b_ij on P(both i and j wrong). Adversary picks any joint law consistent with these. Decision maker picks a deterministic label-symmetric rule to maximize worst-case accuracy.

Claimed results: (A1) optimal rule is weighted majority; (A2) with sum q_i <= 1 weights can be nonnegative; (A4) with sum q_i <= 1, max-min error = min{best single error, B_graph}, with B_graph = inf_y (sum b_ij y_ij)/g(y), g(y) = min over bipartitions of the larger internal edge weight; (A5) the formula fails when sum q_i > 1 (K=5 example).

## 1. Verdict (short)

- No source read solves the pairwise-constrained max-min aggregation problem. The closest explicit statement is de Oliveira, Ishii and Lin (arXiv 2106.00088, Sec. 7.1.2) [FT], which says that overlapping pairwise knowledge of correlations is outside its analysis and that its results "no longer apply". It gives no binary-state pairwise solution.
- The DRO literature on pairwise-known probability bounds (Hunter/Worsley spanning-tree bounds; Boros et al. 2014; Ramachandra and Natarajan 2022) uses the same technical toolkit: LP duality over pairwise-consistent laws, with spanning-tree or star-tree certificates. But its objects differ: it bounds the union from above using pairwise equalities or lower bounds on intersections. Our constraints are upper bounds on co-errors and the event is a threshold (majority) event. The direction is opposite, so I do not think it is the same problem (my reasoning, not a source's claim).
- Robust-aggregation economics (Arieli et al. 2018/2025, Guo et al. 2024, Pan et al. 2023, Neyman–Roughgarden, Babichenko–Garber, Levy–Razin 2022) either treats marginals-only ambiguity or uses different constraints or criteria. None uses pairwise co-error upper bounds with a symmetric majority rule.
- Verdict: to my knowledge, A1/A2/A4/A5 are not in the sources read. Confidence is moderate, because several key items were read only at abstract level (see Sec. 5). The result must cite de Oliveira et al. Sec. 7.1.2 as the explicit open problem, Arieli et al. 2025 as the nearest robust-aggregation result, and the Hunter/Worsley/Boros tree-bound line as the closest technique.

## 2. Source table

| # | Source | Status | Setting | Main result | Overlap with A1/A2/A4/A5 |
|---|---|---|---|---|---|
| 1 | de Oliveira, Ishii, Lin, "Robust Aggregation of Correlated Information", arXiv:2106.00088v3 (econ.TH; metadata shows a 2026-10-04 version; no journal listed) | [FT] | Worst case over all correlation structures among information sources; marginals known | 2 states, 2 actions: use only the best single source. Sec. 7.1.2 handles partial correlation knowledge; no binary pairwise solution | Partial (best-single-source logic parallels the min{best single error, ...} term in A4). Pairwise case: none, explicitly left out |
| 2 | Arieli, Babichenko, Talgam-Cohen, Zabarnyi, "A Random Dictator Is All You Need", arXiv:2302.03667 (v2 Nov 2023); AEJ: Micro 17(1), 2025, pp. 66–96 (per IDEAS) | [ABS] (full-text HTML fetch was refused by the proxy with HTTP 429; not retried) | n symmetric agents, binary recommendations, prior and marginals known, adversarial joint law | For large n, under minimax, regret and approximation-ratio criteria, the unique optimal rule is a random dictator; minimal regret via concavification | Low-partial. Same flavor (robust binary aggregation over adversarial correlation) but marginal-only, not pairwise. Must check: in the marginal-only limit, our framework must reproduce this result |
| 3 | Arieli, Babichenko, Smorodinsky, "Robust forecast aggregation", PNAS 115(52):E12135–E12143 (2018), doi 10.1073/pnas.1813934115 | [ABS] | Binary state; two Blackwell-ordered or conditionally i.i.d. experts; regret with scoring rules | Low-regret schemes for two experts; with many i.i.d. experts, no scheme beats (0.5, 0.5) asymptotically | None to low |
| 4 | Guo, Hartline, Huang, Kong, Shah, Yu, "Algorithmic Robust Forecast Aggregation", arXiv:2401.17743 (Jan 2024) | [ABS] | Finite family of information structures; two agents with independent signals | Approximation schemes for minimax regret | Low (no pairwise constraints) |
| 5 | Pan, Chen, Kong, "Robust Decision Aggregation with Second-order Information", arXiv:2311.14094 (Nov 2023) | [ABS] | Two experts, binary recommendations, conditional independence | Minimax regret; second-order info helps only in specific cases | None |
| 6 | Neyman, Roughgarden, "Are You Smarter Than a Random Expert? The Robust Aggregation of Substitutable Signals", arXiv:2111.03153 (v2 Feb 2022) | [ABS] | Adversarial information structures under a projective-substitutes condition | Averaging beats a random expert; averaging plus extremizing with a known prior | None (no majority or pairwise structure) |
| 7 | Babichenko, Garber, "Learning of Optimal Forecast Aggregation in Partial Evidence Environments", arXiv:1802.07107 (2018); later version Math. Oper. Res. 46(2) (2021), 628–641 (per IDEAS) | [ABS] | Conditionally independent partial evidence; learning from repeated outcomes | Polynomial-time learnability characterization | None to low |
| 8 | Levy, Razin, "Combining forecasts in the presence of ambiguity over correlation structures", JET 199 (2022), doi 10.1016/j.jet.2020.105075 | [ABS] (IDEAS abstract; ScienceDirect not read) | Ambiguity over correlation; correlation bounded via pointwise mutual information; Naive-Bayes benchmark | Predictions set characterized by two parameters; higher correlation bound can make behaviour more conservative | Low-partial in spirit (ambiguity over a set of correlation structures defined by a bound). Different constraint type; no majority rule or pairwise co-error |
| 9 | Balsubramani, Freund, "Optimally Combining Classifiers Using Unlabeled Data", arXiv:1503.01811 (2015) | [ABS] | Binary classifiers; transductive minimax game over unlabeled data; label information as game constraints | Minimax solution indicates when a weighted combination beats the best single classifier | Low-partial (worst-case weighted combination of classifiers). Adversary acts on labels, not joint error laws; no pairwise co-error constraint seen |
| 10 | Ellis, "Condorcet meets Ellsberg", Theoretical Economics 11(3) (2016) | [FT] (PDF via redirect; summarized by fetch tool) | Common-value voters, conditionally independent signals; maxmin over prior of the state only | Plurality voting does not aggregate under ambiguous prior; equilibrium correct-win probability at most 1/2 in some state | Low (ambiguity is over the prior, not correlation; majority is shown to fail, not optimized) |
| 11 | Kaniovski, Zaigraev, "Optimal jury design for homogeneous juries with correlated votes" (WIFO PDF, 2009; later venue not verified) | [FT] | Homogeneous jury; known common pairwise correlation c via Bahadur representation | "Correlation-robust quota" alpha = p(n-1)+1 with floor/ceil choice; robustness = insensitivity to known c, not minimax over c | Low-partial: pairwise correlation constraints in a jury; no worst-case over dependence, weights, or pairwise upper bounds |
| 12 | Boland, Proschan, Tong, "Modelling dependence in simple and indirect majority systems", J. Applied Probability 26(1):81–88 (1989) | [ABS] | Majority systems with dependent components/voters | Monotonicity of reliability functions | Low (cite as majority-with-dependence background; no worst-case result seen) |
| 13 | Kuncheva, Whitaker, "Measures of diversity in classifier ensembles and their relationship with the ensemble accuracy", Machine Learning 51 (2003) | [FT] | Classifier ensembles; empirical study | Double-fault DF = P(both wrong) is the pairwise co-error statistic; diversity correlates weakly with majority-vote gain; no worst-case theorem | Partial in measure only: double-fault equals our b_ij. No worst-case bound. Cite as the ML-side name for our co-error |
| 14 | Ramachandra, Natarajan, "Tight Probability Bounds with Pairwise Independence", arXiv:2006.00516 (v6, 2022) | [FT] | Pairwise independence (equalities P(A_i∩A_j) = p_i p_j) | Closed-form tight union bound min(sum p_i - p_n sum_{i<n} p_i, 1), attained by a star tree; k>=2 bounds via LP duality; threshold events with identical marginals | Low-partial. Same family (worst-case joint law with pairwise info, closed form, tree/star certificate). Equalities, not upper bounds; no majority rule; no non-identical threshold results |
| 15 | Boros, Scozzari, Tardella, Veneziani, "Polynomially Computable Bounds for the Probability of the Union of Events", Math. Oper. Res. (2014) | [NA] for the paper (not fetched). Slide deck by Boros, "Polynomially computable sharp probability bounds" (rutcor.rutgers.edu/Boros.pdf) [FT] | Union probability with known single and pairwise probabilities | Sharp pairwise union upper bound equals best Hunter–Worsley spanning-tree bound (m=2); vertices of the feasible weight region correspond to spanning trees; LP over weight polytopes; inequality version uses lower bounds on intersections | Partial in technique (spanning-tree certificate, LP duality). Different object (union, not threshold); direction differs (union upper bound uses lower bounds on intersections; ours uses upper bounds on co-errors). Highest-priority check |
| 16 | Yang, Alajaji, Takahara, "A Short Survey on Bounding the Union Probability using Partial Information", arXiv:1710.07576 (2017) | [ABS] | Union bounds with individual/pairwise probabilities | Survey; Gallot–Kounias observations | Background only; Hunter/Worsley/graph content not confirmed from abstract |
| 17 | Boros, Lee, "Boole's probability bounding problem, linear programming aggregations, and nonnegative quadratic pseudo-Boolean functions", arXiv:2110.10672 (v4, 2025) | [ABS] | Hailperin LP for Boole's problem; subcones of positive orthant | Hierarchy of aggregations, polyhedral descriptions | Low. No pairwise upper bounds or graph certificates seen in abstract |
| 18 | Miao, Natarajan, Teo, "A Persistency Model and Its Applications in Choice Modeling", Management Science 55(3) (2009), 453–469 | [ABS] (Optimization Online page) | Random objective coefficients; distribution known only via marginals or moments | Persistency probabilities via concave maximization | None |
| 19 | Deza, Laurent, Geometry of Cuts and Metrics (1997); Barahona, Mahjoub (1986) on cut-polytope odd-cycle inequalities | [NA] (book not fetched; cited from background knowledge, verify) | Cut / correlation (Boolean quadric) polytope | Odd-cycle (triangle-type) inequalities describe the cut polytope relaxations | Relevant structurally: pairwise-consistent laws of binary variables form a correlation (Boolean quadric) polytope; odd cycles define its facets. Cited via lecture notes below |
| 19b | Lecture notes "Lectures on 0/1-polytopes", arXiv:math/9909177 | [FT] | Correlation (Boolean quadric) and cut polytopes | Affine equivalence COR(n-1) ~ CUT(n); triangle-type inequalities; odd-cycle inequalities not named; Boole/Pitowsky not covered | Structural background only |
| 20 | Social choice search ("maxmin voting rule correlated voters", "robust jury theorem ambiguity") | Searched twice; no paper solving maxmin voting under correlation constraints found | n/a | Unverified candidates to check: Ladha (1992, 1993) on correlated voters; Berend and Paroush (1998); Dietrich and Spiekermann (2013) | None found |

## 3. Detailed notes on the key items

### 3.1 de Oliveira, Ishii, Lin (arXiv:2106.00088) [FT]
- Section headings (from the full-text fetch): 1 Introduction; 2 Model; 3 Preliminaries (3.1 Payoff polyhedron; 3.2 Blackwell order and supremum; 3.3 Nature's minmax problem); 4 Binary state environment (4.1 binary-state binary-action; 4.2 separable problems; 4.3 general decision problems and decompositions); 5 General-state problems; 6 Behavioral implications (6.1 ignoring information; 6.2 aggregating action recommendations and unanimity); 7 Discussion (7.1 Knowledge of correlation: 7.1.1 Common origin; 7.1.2 Partial knowledge of correlations; 7.2 Ambiguity about marginals; 7.3 Aggregating experts' opinions). The text ends in Sec. 7.3; the conclusion was not in the fetched content.
- Sec. 7.1.2: two cases. (i) Partition blocks: if the joint is known within each block, each block is treated as one marginal experiment. (ii) Overlapping pairs: the agent knows P_12 and P_23 sharing source 2. The text says "our analysis does not immediately extend to other, more complex situations". Feasible set = two pairwise marginal constraints. No binary-state pairwise result, no characterization of an optimal rule, and no mention of majority, weighted majority or graph conditions. Forward-looking references: Pearl (2009) and Spiegler (2016) on correlation restrictions from causal diagrams.
- Use in the paper: cite this section as the explicit open problem. Also note the 2x2 result ("use only the best single source") as the analogue of our min{best single error, ...} term.

### 3.2 Arieli, Babichenko, Talgam-Cohen, Zabarnyi (arXiv:2302.03667; AEJ: Micro 2025) [ABS]
- Abstract-level: n symmetric agents, binary recommendations about a binary state; prior and marginals known; joint law adversarial. For large n, under minimax, regret and approximation-ratio criteria, the unique optimal rule is a random dictator. Minimal regret for any n via concavification.
- Not verified: whether pairwise constraints appear anywhere (no), whether weighted majority is discussed in the body, the exact minimax objective (accuracy vs. expected loss), and whether small n results match our A1/A4.
- Consistency issue to resolve before submission: if the pairwise bounds are vacuous (b_ij = min(q_i, q_j)), does our max-min problem reproduce the random-dictator result? Our A1 says weighted majority is optimal, so either our objective differs (accuracy against state with label symmetry, vs. their objective) or the large-n dictator result applies only outside our sum q_i <= 1 regime. This is the most important reconciliation to make.

### 3.3 Ramachandra and Natarajan (arXiv:2006.00516) [FT]
- Theorem 2.2: tightest union bound under pairwise independence: min(sum_i p_i - p_n sum_{i<n} p_i, 1) with p sorted so p_n is largest. It is the Hunter–Worsley bound with a star tree centered at node n (under pairwise independence).
- Pairwise info is equalities. Upper bounds on pairwise co-occurrence are not treated (the only inequality variant cited is Boros et al. 2014 with lower bounds).
- Threshold events (at least k of n): Thm 3.1 gives ordered Schmidt–Siegel–Srinivasan, Boros–Prékopa and Chebyshev-type bounds; Thm 4.1 gives a tight bound for identical marginals via an aggregated moment LP. Tree structure appears only for k = 1.
- Overlap: low-partial. The pattern "closed-form worst-case over joint laws with pairwise information" is shared; the threshold-event with non-identical marginals and pairwise upper bounds is not covered.

### 3.4 Boros et al. (Math. Oper. Res. 2014) [NA paper; FT slides]
- Slides (Boros, RUTCOR): for m = 2, the sharp pairwise upper bound on union equals the best Hunter–Worsley bound, UB2(FN) = UB2(Fsub) = UB2(HW); the vertices of the weight region are (w1, w2) = (1, X_T) for spanning trees T; the inequality version uses lower bounds on intersections; m = 3 and cherry trees, chordal graphs mentioned; no cycles or odd cycles named; not a correlation-polytope LP.
- Important: "UB2(FN) != UB2*" (pairwise-based bound is not the LP optimum in general). Our B_graph is an inf over weights with a bipartition function g(y); compare carefully to make sure the two are not the same object after LP duality.
- Must do: read the MOR paper to check whether it covers threshold events or upper bounds on pairwise probabilities.

### 3.5 Levy and Razin (JET 2022) [ABS]
- Correlation is bounded through pointwise mutual information; the predicted set is characterized by a Naive-Bayes reading and a correlation bound; the decision maker faces ambiguity over the set. The abstract does not say whether the optimal combination is a weighted rule, nor whether a maxmin over linear rules is optimal. Closest in spirit to "ambiguity over a bounded set of correlation structures", but it is not a pairwise co-error constraint and does not study majority/threshold aggregation.

### 3.6 Kuncheva and Whitaker (2003) [FT]
- Double-fault DF_{i,k} = N00/(N11+N10+N01+N00), the fraction of cases both classifiers get wrong; team value = average over pairs. This is exactly our b_ij statistic (up to normalization).
- No worst-case bound or minimax result in terms of double-fault values; theoretical limits are for averaged Q (Kuncheva et al. 2000).
- Use in the paper: frames pairwise co-error as the known ML diversity measure, and argues that existing ML practice uses averages, not worst-case over joints.

### 3.7 Kaniovski and Zaigraev (WIFO PDF, 2009) [FT]
- Homogeneous jury with a common pairwise correlation c (Bahadur representation; bound on c for valid probabilities). Correlation-robust quota alpha = p(n-1)+1; if alpha is an integer the rule is immune to c, otherwise compare floor/ceil. Robustness is insensitivity to a known c, not a minimax over a set of c.
- Overlap: low-partial. Shows that pairwise correlation constraints are used in the jury literature, with quota rules and no worst-case weight optimization.

### 3.8 Ellis (Theoretical Economics 2016) [FT]
- Ambiguity over the prior of the state only (maxmin EU); signals conditionally independent; plurality voting. The result is nonaggregation under ambiguous priors. Shows that majority-type aggregation is a central object in the ambiguity-in-voting literature, but the ambiguity object is different.

## 4. Sanity check of A5 (done by hand, no code)

Let N = number of wrong votes among K = 5, E[N] = sum q_i = 1.5, and E[C(N,2)] = sum over the 10 pairs of P(both wrong) <= 10 x 0.05 = 0.5.
Identity: C(n,2) - (n-1) = (n-1)(n-2)/2 for every integer n >= 0. Taking expectations: E[C(N,2)] = E[N] - 1 + E[(N-1)(N-2)/2] = 0.5 + E[(N-1)(N-2)/2].
Since (N-1)(N-2)/2 >= 0 for all integers N >= 0, with equality only at N = 1 and N = 2, the pairwise bound forces (N-1)(N-2) = 0 a.s., so N in {1, 2} a.s. Majority of 5 errs iff N >= 3, so P(error) = 0 for every feasible law. This confirms A5 and gives a clean second-moment (Bonferroni-type) certificate. Conclusion: with sum of error rates = 1.5 and pairwise bounds summing to exactly 0.5, the second-moment bound alone is tight enough. This is a known style of argument (second-moment / Bonferroni), so it can be cited as classical, but the specific use for a majority-rule guarantee is what needs to be positioned as new.

## 5. Gaps and what I could not verify

- Arieli et al. 2025 full text: [ABS] only. The proxy returned HTTP 429 on the arxiv.org/html/2302.03667v2 fetch and asked not to refetch that page. Not retried.
- Boros et al. 2014 MOR paper: [NA]; only the slide deck was read.
- Deza and Laurent (1997) and the odd-cycle facets (Barahona and Mahjoub 1986): [NA]; cited from background knowledge.
- Hunter (1976), Worsley (1982), Kounias (1968): cited only as they appear in Ramachandra–Natarajan and Boros slides; not read.
- Guo et al., Pan et al., Arieli et al. 2018, Neyman–Roughgarden, Babichenko–Garber, Levy–Razin, Balsubramani–Freund, Boland et al.: abstract/metadata only.
- Social choice (item 4): two searches only. Ladha (1992, 1993), Berend and Paroush (1998), Dietrich and Spiekermann (2013) are unverified candidates.
- Not searched: DRO work on threshold/chance constraints with pairwise upper bounds in operations research (e.g., Natarajan et al. on chance constraints with pairwise moments); "Doan and Natarajan" was named in the brief but not located in this session.

## 6. Must-cite list (for the paper)

1. de Oliveira, Ishii, Lin (arXiv:2106.00088), Sec. 7.1.2 (explicit open case).
2. Arieli, Babichenko, Talgam-Cohen, Zabarnyi (AEJ: Micro 2025; arXiv:2302.03667) (nearest robust binary aggregation; marginal-only).
3. Arieli, Babichenko, Smorodinsky (PNAS 2018).
4. Levy and Razin (JET 2022) (ambiguity over correlation structures).
5. Boros, Scozzari, Tardella, Veneziani (Math. Oper. Res. 2014); Hunter (1976); Worsley (1982); Kounias (1968) (pairwise-known union bounds, spanning trees).
6. Ramachandra and Natarajan (arXiv:2006.00516) (tight pairwise bounds, LP duality).
7. Kuncheva and Whitaker (Machine Learning 2003) (double-fault = pairwise co-error).
8. Barahona and Mahjoub (1986); Deza and Laurent (1997) (cut/correlation polytope and odd-cycle facets) for B_graph's structure.
9. Ellis (2016) and Kaniovski–Zaigraev (2009) as jury/ambiguity background; Boland et al. (1989) for majority systems with dependence.
10. Guo et al. (2024), Pan et al. (2023), Neyman–Roughgarden (2022), Babichenko–Garber (2021), Balsubramani–Freund (2015): robust aggregation / minimax background.

## 7. Positioning (suggested wording to test)

"Pairwise-constrained robust aggregation has been left open: de Oliveira et al. (Sec. 7.1.2) explicitly note that their analysis does not extend to overlapping pairwise knowledge of correlations, and the robust-aggregation literature (Arieli et al.; Guo et al.; Pan et al.) treats marginal-only ambiguity or other criteria. We solve the binary, label-symmetric case: the optimal rule is a weighted majority (A1); with sum of error rates at most 1 weights can be taken nonnegative (A2); the max-min error equals min{best single error, graph certificate} (A4); and the characterization fails beyond the sum-of-error regime (A5). The technical certificate is a weighted-graph LP dual in the spirit of Hunter/Worsley spanning-tree bounds (Boros et al. 2014; Ramachandra and Natarajan 2022), adapted to upper bounds on co-errors and threshold events."

Do not claim novelty of the Hunter/Worsley technique, second-moment bounds, or the cut/correlation polytope structure.
