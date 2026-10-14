# Max–min weight choice for weighted votes over Fréchet classes: tools from probability, OR and robust aggregation

**Setting (from the project, not re-derived).** Z_i = 1{classifier i correct}, p_i = P(Z_i = 1), q_i = 1 − p_i. A weighted vote without bias is correct iff f_w(Z) = 1{Σ_i w_i(2Z_i − 1) > 0}, a self-dual threshold function of Z. The robust problem is

max_w V(w), where V(w) = min_{π ∈ Π} E_π f_w(Z).

Π is either the Fréchet class (marginals p only) or a smaller set (marginals plus pairwise p_ij, plus correlation bounds, or plus the "even part" of π observable from unlabelled data).

**Notation.**
- [UNVERIFIED] means the item could not be checked against a primary page in this session.
- [INFERENCE] means my own derivation or reasoning, not a sourced claim.

---

## Q1. Worst-case and best-case probability that a weighted sum of dependent Bernoulli variables with given marginals crosses a threshold: closed forms, duals, complexity

### Takeaway
- **Unweighted count ("at least k of n").** The marginals-only problem has a classical closed form: the Rüger (1978) / Morgenstern (1980) formula. It was re-proved constructively in 2025.
- **Weighted sums.** There is no closed form. The problem is a linear program (LP) over the 2^K atoms (Hailperin), and its dual separation step is a knapsack problem.
  - It is NP-hard in general, but only weakly: for two-point marginals ±s_i with probability ½, Kreinovich & Ferson reduce from Partition.
  - It is pseudo-polynomial via DP-based compact LPs (Padmanabhan, Ahipasaoglu, Ramachandra & Natarajan).
- **K = 3 and K = 5.** The LP has 8 or 32 variables, so computation is trivial. Any contribution has to be structural (closed forms, characterisations), not computational.
- **Covering formula [INFERENCE].** For any rule whose failure is monotone in the failure indicators, the worst-case accuracy given only marginals is 1 − min(1, τ*_q), where τ*_q is the q-weighted fractional covering number of the clutter of minimal failure coalitions. This gives Rüger's formula for majority, and gives "best single classifier is robust-optimal" in one line.

### Cited Findings
**Closed forms for counts (k-of-n) given marginals only**
- **Rüger (1978), Morgenstern (1980): tight upper bound for count tails.** Natarajan & Padmanabhan attribute it to both. With p_(1) ≤ … ≤ p_(n) sorted, the bound is U(k) = min( min_{t∈{0,…,k−1}} Σ_{i=1}^{n−t} p_(i)/(k−t), 1 ). The lower bound is obtained by applying U to the complements 1 − p_i. Reference: B. Rüger, "Das maximale Signifikanzniveau des Tests 'Lehne H0 ab, wenn k unter n gegebenen Tests zur Ablehnung führen'", *Metrika* 25:171–178 (1978). — [Natarajan & Padmanabhan, arXiv 1910.06321](https://arxiv.org/pdf/1910.06321)
- **Salako (Dec 2025), constructive DP proof of the generalised Boole–Fréchet bounds.** — [arXiv 2512.09161](https://arxiv.org/pdf/2512.09161)
  - Theorem 3.1: sup P(at least k of n) = min{Σ_{i=1}^{n−r*} p_i/(k−r*), 1}, with p increasing and r* the largest r ≤ k−1 such that p_{n−r+1} ≥ Σ_{i=1}^{n−r+1} p_i/(k−r+1).
  - Corollary 3.3 gives the lower bound max{1 − Σ_{i=r*+1}^{n} (1−p_i)/(n−k+1−r*), 0}.
  - The DP/LP framework accepts any compound event E plus finitely many extra linear constraints. The general LP has 2^n variables.
  - Weighted or threshold events are **not** treated.
- **Hart & Weiss (1997), sharp bound for the Rüger-test region.** P(R) ≤ min_i nα_i/i, sharp for every set of marginals. *Statistics & Probability Letters* 35(1):43–48. — [abstract](https://math.huji.ac.il/~hart/abs/probs.html)
- **Rüschendorf (1981), "Sharpness of Fréchet-bounds".** *Z. Wahrscheinlichkeitstheorie verw. Gebiete* 57:293–302. Related: "Random variables with maximum sums", *Adv. Appl. Probab.* 14:623–632 (1982); "Fréchet bounds and their applications" (1991); "Developments on Fréchet bounds", IMS LN 28:273–296 (1996). — [Rüschendorf publication list](https://www.stochastik.uni-freiburg.de/emeriti/rueschendorf/publications/monge). Content not read [UNVERIFIED beyond bibliographic data].
- **Majority-vote worst case: classifier-fusion literature.**
  - Bax's Theorem 1: for m voters (m odd) with average error p < ½, the worst-case majority error is 2p(1 − 1/(m+1)), approaching "twice the error rate of selecting a classifier at random". — [Bax, arXiv 2104.11833](https://arxiv.org/pdf/2104.11833)
  - Zaigraev & Kaniovski give the same result in competence form, min M = (n(2p−1)+1)/(n+1). — [preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_minmax.pdf)

**Sums of general (two-point / discrete) marginals: duality and complexity**
- **Kreinovich & Ferson (2006) [IJAR 41:331–342]: NP-hardness for sums.** "Computing the best-possible bounds for arbitrary n is an NP-hard … problem."
  - Reduction from Partition, with x_i = ±s_i each with probability ½.
  - The infimum of Prob(Σx_i < 0) equals 0 iff the Partition instance has a solution. — [tech report](https://www.cs.utep.edu/vladik/2004/tr04-25a.pdf)
- **Padmanabhan, Ahipasaoglu, Ramachandra & Natarajan, "Extremal Probability Bounds in Combinatorial Optimization".** DOI 10.1137/21M1442504 (SIAM; journal name inferred from the DOI, probably SIAM J. Optim.). — [arXiv 2109.01591](https://arxiv.org/pdf/2109.01591)
  - Problem: U(r) = max_{θ∈Θ} P_θ(max_{x∈X} c̃ᵀx ≥ r) and L(r) = min_θ, over joints with the given marginals.
  - For a plain sum with support [0,K]: the tight upper bound is an LP polynomial in n and K (Corollary 2.2).
  - For Bernoulli variables the closed form (2.6) is Rüger's.
  - Upper bound over 0/1 V-polytopes: weakly NP-hard, with a pseudo-polynomial algorithm.
  - Lower bound: strongly NP-hard, even restricted to Bernoulli variables.
  - Method: an exponential primal LP, then its dual, with separation by DP over partial sums (equivalence of separation and optimisation, Grötschel–Lovász–Schrijver).
  - Note: their S(w,1) puts weights on **sum values**, Σ_l w_l P(Σc_i = l), **not** on the variables. Weighted-variable Bernoulli sums Σ w_i Z_i are covered only as general two-point marginals {0, w_i}.
- **Embrechts & Puccetti, "Bounds for functions of dependent risks", *Finance & Stochastics* 10(3):341–352 (2006): dual ("standard") bounds for P(ΣX_i < s) with given marginals.** — [IDEAS](https://ideas.repec.org/a/spr/finsto/v10y2006i3p341-352.html). Abstract not visible.
- **Puccetti & Rüschendorf, "Bounds for joint portfolios of dependent risks", *Statistics & Risk Modeling* (2012).** Surveys and improves bounds via rearrangements, duality, conditional moments and reduction. Covers marginals only, overlapping known marginals, and positive-dependence restrictions. — [De Gruyter](https://www.degruyterbrill.com/document/doi/10.1524/strm.2012.1117/html)
- **Embrechts, Wang & Wang (preprint 2014) review of the state of the art.** — [PDF](https://math.uwaterloo.ca/~wang/papers/2015EWW%28FS%29.pdf)
  - Makarov (1981) solves n = 2. Embrechts–Puccetti (2006) give a dual bound for n ≥ 3 with homogeneous marginals.
  - Complete-mixability results by Wang et al. (2013), Puccetti–Rüschendorf (2013) and Bernard et al. (2014).
  - Embrechts–Puccetti–Rüschendorf (2013) introduce the rearrangement algorithm (RA).
  - Worst-case ES is an upper bound for worst-case VaR. With support in [a,b] the gap is at most b − a.
- **Haus, "Bounding stochastic dependence, complete mixability of matrices, and multidimensional bottleneck assignment problems".** Operations Research Letters (DOI 10.1016/j.orl.2014.11.009; volume not verified). — [arXiv 1407.6475](https://arxiv.org/pdf/1407.6475)
  - Theorem 2: "It is strongly NP-complete to decide whether an integral matrix … is completely mixable"; this remains so for a fixed number of columns d.
  - There is a PTAS for fixed d, and pseudo-polynomial time for a fixed number of rows.
  - The swapping algorithm (RA-type) "is not an exact method unless NP ⊆ ZPP".
- **Joint mixability.**
  - Wang & Wang, "Joint mixability" (preprint July 2015; MOR, per the file name, [UNVERIFIED volume]): a tuple (F_1,…,F_n) is jointly mixable if some coupling has a constant sum. Mixability characterises when the lower bound on P(S ≤ s) and the worst-case VaR are attained. Sufficient conditions are given for monotone densities, unimodal-symmetric and elliptical marginals; discrete marginals are barely treated. — [PDF](https://math.uwaterloo.ca/~wang/papers/2015Wang-Wang-MOR.pdf)
  - Koike, Lin & Wang, "Joint mixability and notions of negative dependence", *MOR* 49(4):2786–2802 (2024). — [IDEAS](https://ideas.repec.org/a/inm/ormoor/v49y2024i4p2786-2802.html)

**Self-dual threshold functions (the outer feasible set)**
- **OEIS A001532:** "Number of NP-equivalence classes of self-dual threshold functions of n or fewer variables; number of majority (i.e., decisive and weighted) games with n players": 1, 1, 2, 3, 7, 21, 135, 2470, 175428, 52980624. So K = 3 has 2 types (dictator, majority), K = 5 has 7 and K = 7 has 135. — [OEIS](https://oeis.org/A001532)
- **Chow (1961).** A threshold function is uniquely determined by its degree-0 and degree-1 Fourier coefficients (Chow parameters). O'Donnell & Servedio, "The Chow parameters problem", *SICOMP* 40(1):165–199 (2011), give a PTAS for reconstruction. — [abstract](https://www.cs.columbia.edu/%7Erocco/papers/stoc08.html)

### Inferences
- **[INFERENCE] LP form and duality of the inner problem.**
  - Primal: V(w) = min_π Σ_z π(z) f_w(z) subject to Σ_z π(z) z_i = p_i.
  - Dual: V(w) = max {λ_0 + Σ λ_i p_i : λ_0 + Σ λ_i z_i ≤ f_w(z) ∀z ∈ {0,1}^K}.
  - This is the value at p of the convex closure f_w^− in Dughmi's sense (see Q5).
  - With pairwise information, the minorant becomes a **quadratic** pseudo-Boolean function λ_0 + Σλ_i z_i + Σλ_ij z_i z_j. The dual multipliers are the natural "robust weights" to compare with w (cf. σ* in Balsubramani–Freund, Q4).
- **[INFERENCE] Fractional-covering formula for monotone rules given marginals only.**
  - Let y = 1 − Z be the failure indicators, and let the rule's failure event F be an up-set in y (true for every WMR with w ≥ 0). Then max_π P(F) can be computed with "E y_i ≤ q_i" in place of equality, because raising coordinates keeps you in F. So the dual variables are ≥ 0, and
  - **worst-case accuracy = 1 − min(1, τ*_q)**, where τ*_q = min{Σ μ_i q_i : μ ≥ 0, μ(S) ≥ 1 for every minimal failure coalition S} = max{Σ_S λ_S : Σ_{S∋i} λ_S ≤ q_i, λ ≥ 0}. This is a fractional packing of failure coalitions under capacities q_i.
  - For a self-dual WMR the minimal failure coalitions are the minimal winning coalitions of the game.
  - For k-of-n this LP reduces to Rüger's formula. Its separation step is a min-knapsack (find a coalition S with w(S) ≥ w(N)/2 and μ(S) < 1): weakly NP-hard, pseudo-polynomial, and an FPTAS exists.
  - **Corollary:** μ(S) ≥ 1 forces Σμ_i ≥ 1, so τ*_q ≥ min_i q_i. Hence every monotone rule has worst-case accuracy ≤ max_i p_i, with equality for the dictator of the best classifier. This is a one-line proof, in this special setting, of de Oliveira–Ishii–Lin's best-source theorem (Q4).
- **[INFERENCE] Hardness transfer to weighted votes.** Take Kreinovich–Ferson's reduction with Z_i ~ Bernoulli(½) and weights s_i. Deciding whether min_π P(Σ s_i(2Z_i − 1) > 0) = 0 is NP-hard: E[Σ s_i(2Z_i−1)] = 0, so P(> 0) = 0 forces the sum to be 0 a.s., which is a partition. So exact evaluation of V(w) for general K is weakly NP-hard even with marginals only. Extending this to p_i > ½ (e.g. by padding) is plausible but unproven here.

### Gaps
- No source found giving a **closed form for the worst-case accuracy of a general weighted threshold rule** given marginals. The covering formula above is my derivation and should be checked against the Hailperin/Prékopa literature; it is likely folklore for monotone Boolean functions.
- Embrechts–Puccetti (2006), Puccetti–Rüschendorf (*J. Appl. Probab.* 2013, "Sharp bounds for sums of dependent risks") and Embrechts–Puccetti–Rüschendorf (*J. Banking & Finance* 2013, "Model uncertainty and VaR aggregation") could not be read (paywall/bot-block). Their exact dual formulas are not reproduced here [UNVERIFIED].
- No source found stating NP-hardness specifically for **Bernoulli marginals with arbitrary p_i and weighted sums**. The Kreinovich–Ferson case p_i = ½ is the nearest.

---

## Q2. Discrete moment problems and Boolean–Bonferroni bounds: k-of-n with pairwise information, and anything for weighted counts or threshold functions

### Takeaway
- **Prékopa's binomial-moment LP** (1988, 1990), with closed forms from Boros–Prékopa (1989), gives sharp bounds for "at least / exactly r of n" from aggregated moments S_1 = Σp_i, S_2 = Σp_ij, … This covers simple majority (and the symmetric case) only.
- **Non-aggregated information** (individual p_i, p_ij) leads to Hailperin's exponential LP. Polynomial relaxations exist (Boros–Scozzari–Tardella–Veneziani 2014; Boros–Lee; Qiu–Ahmed–Dey), but computing the exact optimal interval for the union probability is NP-hard (Kaski–Mannila–Mohapatra, ESA 2026). Feasibility of (p_i, p_ij) is NP-complete (Pitowsky 1991).
- **Tractable special cases:** tree-structured pairwise information gives polynomial LPs for P(Σ Z_i ≥ k) (Natarajan–Padmanabhan); pairwise independence gives closed forms for k = 1 (Ramachandra–Natarajan).
- **Weighted counts and threshold functions.** I found **no** result for weighted threshold events with pairwise information. The closest are Yang–Alajaji–Takahara's "partial weighted information" bounds, which are for the union, not for majorities.

### Cited Findings
- **Prékopa, "Boole–Bonferroni inequalities and linear programming", *Operations Research* 36(1):145–162 (1988).** Sharp bounds on P(at least one) from binomial moments via a "simple" LP algorithm. — [IDEAS](https://ideas.repec.org/a/inm/oropre/v36y1988i1p145-162.html)
- **Prékopa, "Sharp bounds on probabilities using linear programming", *Operations Research* 38(2):227–239 (1990).** Sharp bounds for exactly-r and at-least-r from binomial moments, via dual feasible basis theorems. Also gives inequalities among binomial moments and non-consecutive moments. — [IDEAS](https://ideas.repec.org/a/inm/oropre/v38y1990i2p227-239.html)
- **Boros & Prékopa, "Closed form two-sided bounds for probabilities that at least r and exactly r out of n events occur", *MOR* 14(2):317–342 (1989).** Closed-form basis inverses and dual vectors, giving closed-form bounds for dual-feasible bases. — [INFORMS](https://pubsonline.informs.org/doi/fpi/10.1287/moor.14.2.317). Page numbers are from the project bib; the DOI confirms 14(2):317.
- **Mádi-Nagy & Prékopa, "On multivariate discrete moment problems and their applications to bounding expectations and probabilities", *MOR* 29(2):229–258 (2004).** A structural theorem for dual feasible bases when the objective has nonnegative divided differences, with applications to Boolean functions of events. — [summary](https://academicnewsletter.sufe.edu.cn/info/401340)
- **Qiu, Ahmed & Dey, "Strengthened bounds for the probability of k-out-of-n events" (optimization-online, Sept 2013).**
  - Uses "individual probabilities and all joint probabilities of up to m (m < n) events".
  - Strengthens the Prékopa–Gao polynomial-size LP and modifies the objective for general k. Pure LP; no closed form.
  - Authors as listed: S. Ahmed, S. S. Dey, F. Qiu. Journal version (Discrete Applied Mathematics, 2016?) is [UNVERIFIED].
  - — [optimization-online](https://optimization-online.org/2013/09/4043/)
  - **Correction to the brief:** this paper is by Qiu–Ahmed–Dey, not Yang–Alajaji–Takahara.
- **Yang, Alajaji & Takahara.** "On bounding the union probability using partial weighted information" (arXiv 1506.08331; the file name suggests *Stat. Probab. Lett.* 2016 [UNVERIFIED]). Their survey (arXiv 1710.07576) covers Dawson–Sankoff (optimal for aggregated S_1, S_2), de Caen, Kuai–Alajaji–Takahara, Gallot–Kounias, Hunter, and the algorithmic stepwise/greedy bounds. All are for the **union** only. — [survey](https://ar5iv.arxiv.org/html/1710.07576); [weighted info](https://arxiv.org/pdf/1506.08331)
- **Boros, Scozzari, Tardella & Veneziani, "Polynomially computable bounds for the probability of the union of events", *MOR* 39(4):1311–1329 (2014).** Adds constraints to the dual LPs so they become polynomially solvable; three new classes of bounds that dominate most efficient bounds. — [summary](https://academicnewsletter.sufe.edu.cn/info/400961)
- **Boros & Lee, "Boole's probability bounding problem, linear programming aggregations, and nonnegative quadratic pseudo-Boolean functions" (arXiv 2110.10672, v4 2025).** A hierarchy of aggregations of Hailperin's LP, with polyhedral descriptions of the subcones. — [arXiv](https://arxiv.org/abs/2110.10672)
- **Kaski, Mannila & Mohapatra, "Optimal union probability interval is NP-hard", ESA 2026, LIPIcs 388, 129:1–129:19.**
  - Input: probabilities of intersections of some of the events. Computing the optimal interval for P(∪) is NP-hard.
  - This resolves a gap noted by Pitowsky (1991) and Boros et al. (2014, 2026).
  - — [Dagstuhl](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ESA.2026.129)
- **Natarajan & Padmanabhan, "Tree bounds for sums of Bernoulli random variables: a linear optimization approach", *INFORMS J. Optimization* (DOI 10.1287/ijoo.2019.0038).** — [arXiv 1910.06321](https://arxiv.org/pdf/1910.06321)
  - Tight bounds on P(Σ Z_i ≥ k) given p_i and p_ij on a tree are LPs with O(n²) variables and O(n³) constraints.
  - With a complete graph, "verifying if the set Θ is nonempty is an NP-complete problem" (Pitowsky Thm 3.3).
  - Generalises to weighted sums of count probabilities and to order statistics. **Not** to Σ c_i Z_i with variable weights.
  - Kounias's tree bound for k = 1 is Σp_i − Σ_{(i,j)∈T} p_ij.
- **Ramachandra & Natarajan, "Tight probability bounds with pairwise independence" (arXiv 2006.00516).** — [arXiv](https://arxiv.org/abs/2006.00516)
  - Closed-form tight union bound under pairwise independence. Boole/tight ratio ≤ 4/3.
  - New upper bounds for k ≥ 2, tight when all p_i are equal.
  - MOR publication [UNVERIFIED].
- **Kaniovski, "Aggregation of correlated votes and Condorcet's jury theorem", *Theory and Decision* 69(3):453–468 (2010).** — [IDEAS](https://ideas.repec.org:443/a/kap/theord/v69y2010i3p453-468.html)
- **Zaigraev & Kaniovski, "Bounds on the competence of a homogeneous jury", *Theory and Decision* 72(1):89–112 (2012).** Bib details from the project's .tex. — [preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_minmax.pdf)
  - Simple majority with homogeneous p and exchangeable correlations, via the Bahadur representation and the LP over profile probabilities.
  - min M = (n(2p−1)+1)/(n+1); max M = min{2np/(n+1), 1}.
  - The worst-case second-order correlation is (2p−1)/(2p), independent of n.

### Inferences
- **[INFERENCE] Aggregated moments do not suffice for weighted rules.** The binomial-moment machinery (aggregated S_k) works because "at least r" is a symmetric event. A weighted threshold event is not symmetric, so only the non-aggregated Hailperin LP (or Mádi-Nagy–Prékopa-type multivariate moments) applies. The weighted analogue of S_2 in Yang–Alajaji–Takahara's sense, Σ_{i<j} c_ij p_ij with chosen weights, is a natural "aggregated" information model for weighted votes, but it is unstudied for majority events.
- **[INFERENCE] Size for K = 5.** The Fréchet + pairwise polytope in R^{32} has dimension 32 − 1 − 5 − 10 = 16. Its vertices can be enumerated, so a full case analysis (which vertex is least favourable for each of the 7 game types, up to relabelling and sign) is feasible. A closed form like the K = 3 segment (parametrised by t = P(111)) is unlikely to be compact.

### Gaps
- No work found on **sharp bounds for weighted threshold events (weighted majority) given individual and pairwise probabilities**, even for K = 5.
- Could not read Boros–Prékopa's full text to check whether the project's K = 3 symmetric closed form is literally a special case. The project's own .tex flags this too.

---

## Q3. Bounds with partial dependence information: variance/correlation constraints, Bernoulli compatibility and the cut polytope

### Takeaway
- **Risk aggregation.** Adding a variance (or higher-moment) constraint on the sum to known marginals narrows worst-case VaR bounds substantially. These bounds are computed by modified rearrangement algorithms, not closed forms (Bernard–Rüschendorf–Vanduffel 2017; Bernard–Rüschendorf–Vanduffel–Yao 2017 for credit portfolios).
- **Which (p, R) are possible.** The set of attainable Bernoulli correlation/second-moment vectors is a polytope affinely equivalent to the cut polytope (Huber–Marić; Chaganty–Joe; Embrechts–Hofert–Wang). Membership is NP-complete in general (Pitowsky 1991), but complete facet lists exist for small sizes: CUT_5 has 56 facets, CUT_6 368, CUT_7 116,764.
- **Small K is tractable.** For K = 3 or 5, "marginals + correlation bounds" ambiguity sets are explicit small polytopes.
- **ML analogue.** The C-bound and second-order PAC-Bayes bounds are Chebyshev–Cantelli / Markov relaxations that use only the first two moments of the vote margin. They ignore the Bernoulli marginal structure and are not sharp over the Fréchet class.

### Cited Findings
**Variance and moment constraints**
- **Bernard, Rüschendorf & Vanduffel, "Value-at-Risk bounds with variance constraints", *J. Risk & Insurance* 84(3):923–959 (2017).** DOI 10.1111/jori.12108. Abstract not available on the pages reached [content UNVERIFIED]. — [IDEAS](https://ideas.repec.org/a/bla/jrinsu/v84y2017i3p923-959.html); [Rüschendorf list](https://www.stochastik.uni-freiburg.de/emeriti/rueschendorf/publications/financial)
- **Bernard, Rüschendorf, Vanduffel & Yao, "How robust is the value-at-risk of credit risk portfolios?", *European J. Finance* 23(6):507–534 (2017).**
  - Loan portfolios with known marginal loss distributions, i.e. Bernoulli-type defaults.
  - Modifies the RA "to approximate sharp VaR bounds when … higher-order aggregate moments, such as variance and skewness, are also known".
  - Moment information "significantly tightens" the bounds, but high-confidence VaR stays non-robust.
  - — [IDEAS](https://ideas.repec.org/a/taf/eurjfi/v23y2017i6p507-534.html)
- **Bernard, Rüschendorf, Vanduffel & Wang, "Risk bounds for factor models", *Finance & Stochastics* 21(3):631–659 (2017).** Each risk has a known joint distribution with a common factor, but conditional independence is dropped. "The dependence uncertainty spread is typically reduced substantially." — [IDEAS](https://ideas.repec.org/a/spr/finsto/v21y2017i3d10.1007_s00780-017-0328-4.html)
- **Rüschendorf & coauthors, other partial-information papers.** — [list](https://www.stochastik.uni-freiburg.de/emeriti/rueschendorf/publications/financial)
  - "Risk bounds with additional information on functionals of the risk vector", *Dependence Modeling* 6:102–113 (2018).
  - Rüschendorf & Witting, "VaR bounds in models with partial dependence information on subgroups", *Dependence Modeling* 5(1):59–74 (2017).
  - "Improved Hoeffding–Fréchet bounds and applications to VaR estimates" (Springer 2017).
- **Bartl, Kupper, Lux, Papapantoleon & Eckstein, "Marginal and dependence uncertainty: bounds, optimal transport, and sharpness" (SIAM journal, DOI 10.1137/21M144709X).** OT duality for Fréchet classes with extra information. Improved Fréchet–Hoeffding bounds are pointwise sharp under marginal uncertainty, but not without it. — [arXiv 1709.00641](https://arxiv.org/pdf/1709.00641)
- **Dhara, Das & Natarajan, "Worst-case expected shortfall with univariate and bivariate marginals", *INFORMS J. Computing* 33(1):370–389 (2021).** Known univariate marginals plus a KL ball around bivariate marginals. Efficiently computable by convex optimisation when the bivariate information is tree-structured. — [arXiv 1701.04167](https://arxiv.org/abs/1701.04167); [IDEAS](https://ideas.repec.org/a/inm/orijoc/v33y2021i1p370-389.html)
- **Vandenberghe, Boyd & Comanor, "Generalized Chebyshev bounds via semidefinite programming", *SIAM Review* 49(1):52–64 (2007).** A sharp lower bound on the probability of a set defined by quadratic inequalities, given first and second moments, via SDP. This is the moment-only (no Bernoulli support) analogue. — [Boyd page](https://stanford.edu/~boyd/papers/prob_bnds.html)

**Bernoulli compatibility and the cut polytope**
- **Chaganty & Joe, "Range of correlation matrices for dependent Bernoulli random variables", *Biometrika* 93(1):197–206 (2006).** Necessary and sufficient conditions for compatibility of (p, R). Some R are compatible with no p. Odds ratios are "less constrained" than correlations. — [IDEAS](https://ideas.repec.org/a/oup/biomet/v93y2006i1p197-206.html)
- **Huber & Marić, "Bernoulli correlations and cut polytopes" (arXiv 1706.06182, 2017).** — [arXiv](https://arxiv.org/abs/1706.06182)
  - For symmetric Bernoulli variables, the set R(B_n) of correlation vectors is a polytope with R(B_n) = 1 − 2·CUT(n).
  - Its vertices are the correlation vectors of uniform distributions on cube diagonals.
  - Realisability can be decided by LP. A sampler is given.
- **Embrechts, Hofert & Wang, "Bernoulli and tail-dependence compatibility" (DOI 10.1214/15-AAP1128, *Ann. Appl. Probab.*).** A unit-diagonal matrix is a tail-dependence matrix iff it is a Bernoulli-compatible matrix (E[XXᵀ] for Bernoulli X) times a constant. — [arXiv 1606.08212](https://arxiv.org/abs/1606.08212)
- **Pitowsky, "Correlation polytopes: their geometry and complexity", *Math. Programming* 50:395–414 (1991).** "Deciding membership in a correlation polytope is an NP-complete problem." — [IAOR](https://iaorifors.com/paper/8045)
- **Shyamalkumar & Tao.** TDM/Bernoulli realisation can be written as an LP whose constraint structure relates to max-cut. Symmetry and sparsity give polynomial tests for parametric families. — [arXiv 1901.02157](https://arxiv.org/abs/1901.02157v2)
- **Facet counts of the cut polytope (SMAPO).** CUT(5): 56 facets in 2 classes. CUT(6): 368 in 3. CUT(7): 116,764 in 11. CUT(8): 217,093,472 (conjectured complete). — [SMAPO](https://comopt.ifi.uni-heidelberg.de/software/SMAPO/cut/cut.html)
- **Boros & Hammer, "Cut-polytopes, Boolean quadric polytopes and nonnegative quadratic pseudo-Boolean functions", *MOR* 18(1):245ff (1993).** Title and DOI seen in search; authors and content [UNVERIFIED]. This is the duality between the Boolean quadric (= correlation) polytope and nonnegative quadratic pseudo-Boolean functions, which is exactly the dual of the "marginals + pairwise" problem. — [INFORMS](https://pubsonline.informs.org/doi/10.1287/moor.18.1.245)
- **Fontana & Semeraro.**
  - "Geometrical representation and dependence structure of three-dimensional Bernoulli distributions" (arXiv 2512.06786, Dec 2025): all extremal points of the 3-dimensional Bernoulli class with common mean p, in closed form. There are 6 for p ≤ 1/3 and 9 for 1/3 < p ≤ ½, each with at most 4 nonzero probabilities. Pairwise correlations attain their bounds at extremal points. — [arXiv](https://arxiv.org/html/2512.06786v1)
  - "High dimensional Bernoulli distributions: algebraic representation and applications" (*Bernoulli* 30(1), 2024, per the URL): polynomial generators of the class. — [arXiv 2205.12744](https://arxiv.org/abs/2205.12744)
  - Fontana, Luciano & Semeraro, "Model risk in credit risk": VaR bounds attained on ray densities of exchangeable classes E_d(p) and E_d(p, ρ). — [ar5iv](https://ar5iv.arxiv.org/html/1906.06164)

**ML second-order bounds (moment relaxations)**
- **Germain, Lacasse, Laviolette, Marchand & Roy, "Risk bounds for the majority vote: from a PAC-Bayesian analysis to a learning algorithm", *JMLR* 16:787–860 (2015).** The C-bound depends on "the average quality of the voters and their average disagreement". MinCq minimises it by QP. — [JMLR](https://jmlr.org/papers/v16/germain15a.html)
  - [UNVERIFIED formula from memory]: C = 1 − (E M)²/E[M²] for the margin M. This is the one-sided Chebyshev (Cantelli) bound.
- **Masegosa, Lorenzen, Igel & Seldin, "Second order PAC-Bayesian bounds for the weighted majority vote" (NeurIPS 2020).** — [ar5iv](https://ar5iv.labs.arxiv.org/html/2007.13532)
  - L(MV_ρ) ≤ 2E_ρ[L(h)] (first-order); L(MV_ρ) ≤ 4E_{ρ²}[L(h,h′)] (tandem loss).
  - Binary case: L(MV_ρ) ≤ 4E_ρ[L(h)] − 2E_{ρ²}[D(h,h′)], where the disagreements "can be estimated from unlabeled data".

### Inferences
- **[INFERENCE] Translation to votes.** Knowing all pairwise correlations fixes Var(Σ w_i(2Z_i−1)) for **every** w. Knowing only an aggregate (average correlation, or the variance of the margin for a fixed w) is the analogue of Bernard–Rüschendorf–Vanduffel's variance constraint. "Marginals + one quadratic constraint" is an LP with one extra row, so its sharp value is a 2-multiplier dual (affine minorant + one quadratic term).
- **[INFERENCE] Sharpness gap of the C-bound.** The C-bound / Cantelli is the sharp worst case over **all** distributions of the margin with given mean and second moment. The project's Fréchet + pairwise version adds Bernoulli-support and marginal constraints, so its sharp value is at least as large and generally strictly better. Quantifying this gap for K = 3 and 5 is a clean, publishable comparison.
- **[INFERENCE] Correlation-box ambiguity sets.** For K = 3 the feasible (p_i, p_ij) set is COR_3 (affinely CUT_4, described by triangle-type/Boole inequalities). For K = 5 it is COR_5 ≅ CUT_6 with 368 facets; this covariance-mapping identity is from my recollection of Deza–Laurent [UNVERIFIED]. An ambiguity set {|ρ_ij − ρ̂_ij| ≤ ε} ∩ COR is therefore an explicit polytope for K ≤ 6.

### Gaps
- Abstract and main theorem of Bernard–Rüschendorf–Vanduffel (JRI 2017) not seen. Whether their method gives sharp (not just RA-approximate) bounds for **discrete/Bernoulli** marginals is unverified.
- No work found that gives worst-case accuracy of **weighted votes** under correlation-ball or variance constraints, or that compares the C-bound with the sharp Fréchet + second-moment bound.

---

## Q4. Robust aggregation in economics and CS theory: how they formalise unknown information structures, what max–min solutions they find, and whether binary votes with known accuracies are covered

### Takeaway
- **Known marginals only: solved.** Following the single best source is robust-optimal for binary-state, binary-action problems (de Oliveira–Ishii–Lin). With symmetric agents, a random dictator is optimal for large n under maximin, regret and approximation-ratio criteria (Arieli–Babichenko–Talgam-Cohen–Zabarnyi). Simple majority is "far from optimal". These are exactly the "votes with known accuracies, adversarial correlation" models; the answer is "do not aggregate".
- **Non-trivial aggregation needs structural restrictions.** Examples: conditional independence, Blackwell ordering, projective substitutes, or extra observables such as second-order predictions (Arieli–Babichenko–Smorodinsky; Neyman–Roughgarden; Pan–Chen–Kong; Guo et al.; Chen–Peng–Tang 2026).
- **"Even-part" observables.** The closest ML analogue of "marginals + even-part (unlabelled) information" is Balsubramani–Freund. With the vote distribution known and accuracy lower bounds, the minimax-optimal predictor is a clipped weighted vote with weights from a convex dual.

### Cited Findings
- **de Oliveira, Ishii & Lin, "Robust aggregation of correlated information" (arXiv 2106.00088 v2, Sept 2024; earlier "Robust Merging of Information").** — [arXiv HTML](https://arxiv.org/html/2106.00088v2); [PDF](https://arxiv.org/pdf/2106.00088)
  - Model: m sources with known marginal experiments P_j: Θ → ΔY_j. Joint set J = {P : Σ_{y_−j} P(y|θ) = P_j(y_j|θ)}. Value V = max_σ min_{P∈J} Σ_θ Σ_y P(y|θ) u(θ, σ(y)), with σ: Y → Δ(A) allowed to randomise.
  - **Theorem 1 (|Θ| = |A| = 2):** a best-source strategy is robustly optimal, and it is unique if the best source is unique.
  - Proof idea: the Blackwell supremum's feasible set is the convex hull of the individual feasible sets (Lemma 2), and linearity puts the maximum at an extreme point.
  - Theorem 2: binary state with general actions decomposes into n−1 binary comparisons, using at most n−1 sources.
  - Theorem 3: general states need at most N sources.
  - The minimax theorem is invoked: the objective is bilinear in σ and P, and both sets are convex and compact.
  - Section 6.1 contains **no** analysis of extra correlation knowledge. It is an example in which two sources beat one "regardless of the correlations".
- **Arieli, Babichenko, Talgam-Cohen & Zabarnyi, "A random dictator is all you need", *AEJ: Microeconomics* 17(1):66–96 (2025).** Also circulated as "Universally robust information aggregation for binary decisions", arXiv 2302.03667. — [IDEAS](https://ideas.repec.org/a/aea/aejmic/v17y2025i1p66-96.html); [ar5iv](https://ar5iv.labs.arxiv.org/html/2302.03667)
  - Model: symmetric agents with binary signals. The DM knows the prior μ and each agent's marginal (posteriors p_1 < ½ < p_2), but the correlation is adversarial. The DM sees the fraction of H reports.
  - For many agents, random dictator is uniquely optimal under maximin, regret and approximation ratio.
  - Minimal regret equals the maximal distance between P* and its concavification (Theorem 5.1).
  - Random dictator is not optimal for n = 2 (Example 5.2).
- **Arieli, Babichenko & Smorodinsky, "Robust forecast aggregation", *PNAS* 115(52):12135–12143 (2018).** Regret against the omniscient Bayesian aggregator, binary state. Low-regret schemes exist for two Blackwell-ordered or conditionally i.i.d. experts. With many conditionally i.i.d. experts, nothing beats the uninformative (0.5, 0.5) forecast asymptotically. — [IDEAS](https://ideas.repec.org/a/nas/journl/v115y2018pe12135-e12143.html)
- **Pan, Chen & Kong, "Robust decision aggregation with second-order information" (arXiv 2311.14094; WWW 2024, DOI 10.1145/3589334.3645384).** — [arXiv](https://arxiv.org/html/2311.14094v1)
  - Model: two experts with binary recommendations, plus each expert's prediction of the other's recommendation. Minimax regret.
  - General structures: regret ≥ 0.5, and "follow the first expert" is optimal.
  - Conditional independence: a deterministic threshold aggregator achieves 1/3 using the second-order information.
  - Homogeneous experts: 3 − 2√2.
  - Extension to more experts is left as future work.
- **Neyman & Roughgarden.**
  - "From proper scoring rules to max-min optimal forecast aggregation", *Operations Research* 71(6):2175–2195 (2023): quasi-arithmetic pooling maximises worst-case (over outcomes) profit relative to paying experts by their weights. Linear pooling for the quadratic rule, logarithmic for log. — [IDEAS](https://ideas.repec.org/a/inm/oropre/v71y2023i6p2175-2195.html); [ar5iv](https://ar5iv.labs.arxiv.org/html/2102.07081)
  - "Are you smarter than a random expert? The robust aggregation of substitutable signals" (arXiv 2111.03153): without restrictions no nontrivial guarantee is possible. Under "projective substitutes", averaging beats a random expert, and with the prior known, averaging plus extremising does better still. Venue [UNVERIFIED; likely SODA 2022]. — [arXiv](https://arxiv.org/abs/2111.03153)
- **Guo, Hartline, Huang, Kong, Shah & Yu, "Algorithmic robust forecast aggregation" (arXiv 2401.17743).** Efficient approximation schemes for minimax-regret aggregators when the family of information structures is finite, and for the two-agent conditionally independent setting. — [arXiv](https://arxiv.org/abs/2401.17743)
- **Chen, Peng & Tang, "Prior-agnostic robust forecast aggregation" (arXiv 2604.24517, Apr/May 2026).** A closed-form log-odds (logit-linear) aggregator.
  - Regret < 0.0226 under conditional independence with a known binary state.
  - With known expert marginals, a generalised log-odds rule gets 0.0228 against a lower bound of 0.0225.
  - — [arXiv](https://arxiv.org/abs/2604.24517)
- **Babichenko & Garber, "Learning optimal forecast aggregation in partial evidence environments", *MOR* 46(2):628–641 (2021).** — [IDEAS](https://ideas.repec.org/a/inm/ormoor/v46y2021i2p628-641.html)
- **Levy & Razin, "Combining forecasts in the presence of ambiguity over correlation structures", *JET* 199:105075 (2022).** Bounded correlation measured by pointwise mutual information. The prediction set is determined by the naive-Bayes reading plus the correlation bound. — [IDEAS](https://ideas.repec.org/a/eee/jetheo/v199y2022ics0022053118303144.html)
- **Kaniovski & Zaigraev, "Optimal jury design for homogeneous juries with correlated votes", *Theory and Decision* 71(4):439–459 (2011).** A "correlation-robust" quota α = p(n−1)+1 makes accuracy independent of the (exchangeable) correlation coefficient. This is not max–min. — [IDEAS](https://ideas.repec.org/a/kap/theord/v71y2011i4p439-459.html); [preprint](https://Serguei.Kaniovski.wifo.ac.at/fileadmin/pdf/condorcet_robust.pdf)
- **Balsubramani & Freund, "Optimally combining classifiers using unlabeled data" (COLT 2015, arXiv 1503.01811).** — [arXiv PDF](https://arxiv.org/pdf/1503.01811)
  - Game: V = max_g min_{z∈[−1,1]^n, (1/n)Fz ≥ b} (1/n)zᵀg.
  - Theorem 2: V = −γ(σ*), with γ(σ) = (1/n)Σ_j[|x_jᵀσ| − 1]_+ − bᵀσ minimised over σ ≥ 0. The optimal g_j* = clip(x_jᵀσ*): a weighted vote, randomised inside the "zero box".
  - Lemma 7 gives the dual of the inner problem.
  - Without F, the solution is the single best classifier. With F, gains can be p-fold.
  - Extension: "Optimal binary classifier aggregation for general losses", NeurIPS 2016 (sigmoid-like link of a weighted vote). — [NeurIPS](https://proceedings.neurips.cc/paper/2016/hash/eaa52f3366768bca401dca9ea5b181dd-Abstract.html)
- **Unlabelled-agreement identifiability (relevant to "even-part" information).**
  - Parisi, Strino, Nadler & Kluger (PNAS 2014): under conditional independence, the off-diagonal covariance of predictions is rank-one, with entries proportional to balanced accuracies. — [arXiv 1303.3257](https://arxiv.org/abs/1303.3257)
  - Jaffe, Nadler & Kluger (AISTATS 2015). — [PMLR](https://proceedings.mlr.press/v38/jaffe15.html)
  - Platanios, Blum & Mitchell (UAI 2014), agreement rates. — [PDF](https://www.auai.org/uai2014/proceedings/individuals/313.pdf)
  - Corrada-Emmanuel (arXiv 2312.05392), algebraic "agreement equations" for pair-correlated binary classifiers. — [arXiv](https://arxiv.org/abs/2312.05392)
- **Nitzan & Paroush, "Optimal decision rules in uncertain dichotomous choice situations", *International Economic Review* 23(2):289–297 (1982).** Log-odds weights under independence. — [IDEAS](https://ideas.repec.org/a/ier/iecrev/v23y1982i2p289-97.html)
  - Berend & Kontorovich restate the optimal rule as w_i = log(p_i/(1−p_i)), with error ≤ exp(−Φ/2) where Φ = Σ(p_i − ½)w_i. They call independence "truly restrictive". — [ar5iv](https://ar5iv.arxiv.org/html/1312.0451)

### Inferences
- **[INFERENCE] How the information models compare.**
  - de Oliveira et al.: class-conditional marginal experiments known.
  - Arieli et al.: posterior-type marginals plus prior known.
  - This project: unconditional accuracies p_i = P(Z_i = 1) known.
  - For neutral (self-dual) rules, correctness depends only on Z, so the project's model is a coarsening. The best-source conclusion carries over; it is consistent with the covering-formula corollary in Q1 and with the companion note's numerical check.
- **[INFERENCE] "Even-part" information as a population Balsubramani–Freund game.**
  - Unlabelled data reveal the law of the vote vector X = Y·(2Z − 1). If Y is symmetric and independent of the correctness pattern (e.g. a balanced label with a label-symmetric noise model), this law equals the symmetrisation ½(π(z) + π(1−z)), the even part of π.
  - Adding accuracy equalities gives Balsubramani–Freund's game in population form, with equality constraints. The dual σ becomes sign-free, which permits the negative-weight rules seen in the K = 3 experiments in the companion note.
  - For general K, the even part contains 2^{K−1} − 1 free numbers, versus K(K−1)/2 pairwise agreements. So "even part" ⊋ "pairwise" once K ≥ 4.
- **[INFERENCE] Shape of a new contribution.** Every positive aggregation result above comes from **restricting the dependence**. A clean contribution is to restrict through **observable second-order statistics** (pairwise success probabilities, or agreement rates from unlabelled data) rather than structural assumptions (CI, Blackwell order). That middle ground is unoccupied for K ≥ 3 voters.

### Gaps
- None of the robust-aggregation papers treats **K ≥ 3 binary votes with known pairwise joint accuracies**. Pan–Chen–Kong treat two experts. Arieli et al. treat symmetric agents with marginals only.
- de Oliveira–Ishii–Lin's Section 6 claims some extensions to "additional knowledge about correlations", but the text read does not specify them.
- Venues not verified: Neyman–Roughgarden 2111.03153 and Guo et al. 2401.17743.

---

## Q5. Distributionally robust optimisation over Fréchet classes and multi-marginal OT: max–min of threshold-event probabilities, and which minimax theorems apply

### Takeaway
- **The inner problem is a convex/concave closure.** For a fixed rule, the inner minimum over the Fréchet class is the convex closure f^−(p) of the rule's correctness function, a standard object in submodular optimisation and correlation-gap analysis. Evaluating such closures is NP-hard in general (concave closure of cut functions).
- **Tractability needs structure:** supermodularity, tree or non-overlapping marginals, or small support.
- **Minimax over randomised rules holds by von Neumann/Sion:** the payoff is bilinear on finite-dimensional compact convex sets. So the robust randomised rule is a Bayes rule against a least-favourable joint distribution π*. Deterministic weighted votes can be strictly worse; no general result on that gap was found.

### Cited Findings
- **Dughmi, "Submodular functions: extensions, distributions, and algorithms. A survey" (arXiv 0912.0322).** — [arXiv PDF](https://arxiv.org/pdf/0912.0322)
  - Convex closure f^−(x): the expectation of f under the marginal-x distribution that minimises E f. Concave closure f^+ uses the maximising distribution.
  - "When f is submodular, f^− = L_f" (the Lovász extension), so it is efficiently computable.
  - Theorem 3.11: "It is NP-hard to evaluate f^+(x) for an arbitrary submodular f", even for graph cut functions.
- **Agrawal, Ding, Saberi & Ye, "Correlation robust stochastic optimization" / "Price of correlations in stochastic optimization", *Operations Research* (DOI 10.1287/opre.1110.1011).** — [ar5iv](https://ar5iv.arxiv.org/html/0902.1792)
  - Marginals known, correlation unknown; minimise the worst-case expected cost.
  - The correlation gap is the "ratio κ between the expected cost of the worst case distribution and that of the independent distribution".
  - κ ≤ e/(e−1) for submodular costs, and can be Ω(2^n) for some supermodular functions.
  - For supermodular costs the worst case is analytic, giving polynomial algorithms.
- **Doan & Natarajan, "On the complexity of nonoverlapping multivariate marginal bounds for probabilistic combinatorial optimization problems", *Operations Research* 60(1):138–149 (2012).** Tight bounds given non-overlapping multivariate marginals, polynomial in some discrete cases. The bound is exactly attained. — [IDEAS](https://ideas.repec.org/a/inm/oropre/v60y2012i1p138-149.html)
- **Natarajan, Song & Teo, "Persistency model and its applications in choice modeling", *Management Science* 55(3):453–469 (2009).** Marginal-distribution model for discrete optimisation. Content not read. — [IDEAS](https://ideas.repec.org/a/inm/ormnsc/v55y2009i3p453-469.html)
- **Gao & Kleywegt, "Distributionally robust stochastic optimization with dependence structure" (arXiv 1701.04200).** Wasserstein plus second-moment constraints, and rank-dependence/copula Wasserstein sets, with dual reformulations. Not a fixed-marginal formulation. — [arXiv](https://arxiv.org/abs/1701.04200)
- **Chen, Sim & Xiong, "Robust stochastic optimization made easy with RSOME", *Management Science* 66(8):3329–3339 (2020).** Event-wise ambiguity sets; a modelling tool. — [IDEAS](https://ideas.repec.org/a/inm/ormnsc/v66y2020i8p3329-3339.html)
- **Vovk & Wang, "Combining p-values via averaging" (arXiv 1212.4966; *Biometrika* 107(4):791–808).** Twice the arithmetic mean is valid under arbitrary dependence: "an old result by Rüschendorf", and "the factor of 2 cannot be improved". Merging functions are characterised via worst-case quantiles over all dependence structures (Lemma 1). Weighted averages are only "barely" treated. — [ar5iv](https://ar5iv.arxiv.org/html/1212.4966)
  - Follow-up: Vovk, Wang & Wang, "Admissible ways of merging p-values under arbitrary dependence", *Ann. Statist.* 50(1) (2022). — [arXiv 2007.14208](https://arxiv.org/pdf/2007.14208)
- **de Oliveira–Ishii–Lin** explicitly use the minimax theorem on the bilinear σ–P game (see Q4). — [arXiv](https://arxiv.org/pdf/2106.00088)

### Inferences
- **[INFERENCE] Minimax over randomised rules.**
  - Let G be the set of randomised neutral rules g: {±1}^K → [0,1] with g(−x) = 1 − g(x); this is a polytope. Let Π be the ambiguity polytope. Then max_{g∈G} min_{π∈Π} A(g, π) = min_{π∈Π} max_{g∈G} A(g, π), by von Neumann, since A is bilinear.
  - Hence the robust randomised rule is Bayes-optimal for some least-favourable π*. Under any π*, the Bayes rule is a likelihood-ratio rule on the vote profile.
  - If π* turned out to be conditionally independent, the rule would be a log-odds WMR. This gives a concrete route to "robust log-odds weights": characterise π*.
  - For deterministic WMRs the max–min is over a finite set (7 game types up to symmetry for K = 5). The deterministic/randomised gap is an open quantitative question. The companion note found it numerically zero for K = 3 with pairwise information.
- **[INFERENCE] Not a supermodular problem.** Majority-type correctness functions are neither submodular nor supermodular, so the Agrawal et al. analytic worst-case (supermodular) does not apply directly. The correlation gap of a WMR (worst case versus independence) is a natural new quantity. Its asymptotics in K for log-odds weights are unstudied.

### Gaps
- No DRO paper found that maximises over **decision rules** the worst-case probability of a threshold event over a Fréchet class with pairwise constraints. Existing DRO work evaluates fixed objectives or optimises continuous decisions.
- Natarajan–Song–Teo and Chen–Sim–Xiong were not read beyond bibliographic data.

---

## Tools we can use and what is open

### Takeaway
Most building blocks exist: Hailperin/Prékopa LP duality, the Rüger closed form, tree-structured LPs, the cut-polytope geometry, Balsubramani–Freund's dual, and von Neumann minimax. With **marginals only** the max–min is settled: the best single classifier wins. The open, plausibly solvable space is max–min **weight design under observable second-order information** (pairwise accuracies or agreements, correlation balls, or the even part), for K = 5 and for general K with sparse or tree-structured dependence.

### Cited Findings (directly applicable results)
1. **Hailperin LP and its dual.** Worst-case accuracy of any rule is an LP over 2^K atoms. Its dual is the best affine (marginals) or quadratic (pairwise) pseudo-Boolean minorant. K = 3 and 5 are trivial to compute. — [Natarajan–Padmanabhan](https://arxiv.org/pdf/1910.06321); [Salako](https://arxiv.org/pdf/2512.09161); [Boros–Lee](https://arxiv.org/abs/2110.10672)
2. **Rüger/Morgenstern closed form** for simple majority given marginals; constructive DP proof by Salako. — [Natarajan–Padmanabhan](https://arxiv.org/pdf/1910.06321); [Salako](https://arxiv.org/pdf/2512.09161)
3. **Prékopa / Boros–Prékopa** sharp bounds from aggregated binomial moments: symmetric K, simple majority, closed forms. — [Prékopa 1990](https://ideas.repec.org/a/inm/oropre/v38y1990i2p227-239.html); [Boros–Prékopa 1989](https://pubsonline.informs.org/doi/fpi/10.1287/moor.14.2.317)
4. **Tree-structured pairwise information gives a polynomial LP** for count thresholds. — [Natarajan–Padmanabhan](https://arxiv.org/pdf/1910.06321)
5. **Pseudo-polynomial compact LPs via DP separation for weighted sums**, and NP-hardness of exact evaluation. — [Padmanabhan et al.](https://arxiv.org/pdf/2109.01591); [Kreinovich–Ferson](https://www.cs.utep.edu/vladik/2004/tr04-25a.pdf)
6. **Marginals only: best source is robust-optimal; random dictator for symmetric agents.** — [de Oliveira–Ishii–Lin](https://arxiv.org/pdf/2106.00088); [Arieli et al.](https://ideas.repec.org/a/aea/aejmic/v17y2025i1p66-96.html)
7. **Known vote distribution plus accuracy bounds: clipped weighted vote**, with weights from the convex dual γ(σ). — [Balsubramani–Freund](https://arxiv.org/pdf/1503.01811)
8. **Compatibility geometry.** COR/CUT polytope; NP-complete in general; full facet lists for small n. — [Huber–Marić](https://arxiv.org/abs/1706.06182); [Pitowsky](https://iaorifors.com/paper/8045); [SMAPO](https://comopt.ifi.uni-heidelberg.de/software/SMAPO/cut/cut.html); [Chaganty–Joe](https://ideas.repec.org/a/oup/biomet/v93y2006i1p197-206.html)
9. **Moment-constrained risk bounds.** Variance and skewness constraints via modified RA; factor models. — [BRVY 2017](https://ideas.repec.org/a/taf/eurjfi/v23y2017i6p507-534.html); [BRVW 2017](https://ideas.repec.org/a/spr/finsto/v21y2017i3d10.1007_s00780-017-0328-4.html)
10. **Hardness landmarks for general K with pairwise information.** — [Pitowsky](https://iaorifors.com/paper/8045); [Kaski–Mannila–Mohapatra ESA 2026](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ESA.2026.129); [Dughmi](https://arxiv.org/pdf/0912.0322)

### Inferences (tools derived here, to be checked)
- **[INFERENCE] Covering formula.** For a monotone rule given marginals only, worst-case accuracy = 1 − min(1, τ*_q), where τ*_q is the fractional covering number of the minimal failure (winning) coalitions with costs q_i = 1 − p_i.
  - It gives Rüger's formula for majority, and "≤ max_i p_i" for every monotone rule in one line.
  - Separation is a min-knapsack: pseudo-polynomial, with an FPTAS.
  - Likely folklore, so check before claiming.
- **[INFERENCE] Minimax over randomised rules** reduces robust weight design to finding a least-favourable π*. Robust weights can then be read off π*'s likelihood ratios or the LP dual multipliers.

### Gaps (what appears open, roughly in order of tractability)
1. **K = 5, marginals + full pairwise p_ij: robust-optimal WMR.**
   - Characterise, as a function of (p, p_ij), which of the 7 game types (with signs) is max–min optimal.
   - Determine when it beats the best single classifier.
   - Compare it with log-odds weights.
   - The Fréchet + pairwise polytope has dimension 16. No published treatment found; K = 3 is done in the project.
2. **Dual / structural characterisation for general K.**
   - Robust value = max over quadratic pseudo-Boolean minorants of f_w.
   - Open: whether the optimal w can be read off the dual multipliers (as with σ* in Balsubramani–Freund); whether π* is ever conditionally independent; whether signed (negative) weights are necessary.
3. **Complexity with pairwise information.** Exact evaluation is NP-hard for general K (Pitowsky). Open:
   - tree or sparse dependence graphs for weighted thresholds: a plausible pseudo-polynomial DP extending Natarajan–Padmanabhan;
   - approximation guarantees;
   - hardness of the **outer** max over w (not just the inner evaluation).
4. **Correlation-ball and variance ambiguity sets.**
   - Robust weights as ε grows from 0 (point-identified π) to ∞ (Fréchet class, i.e. dictator): a phase-transition / interpolation result.
   - The analogue of Bernard–Rüschendorf–Vanduffel for votes.
   - The sharp Fréchet + second-moment bound against the C-bound / tandem bound.
5. **Even-part (unlabelled) information with equality accuracy constraints.**
   - Population version of Balsubramani–Freund. Sign-free σ explains negative weights.
   - K ≥ 4 with **only pairwise agreements** (not the full vote law) is unstudied.
6. **Statistical layer.** Robust weights when p and p_ij are estimated: confidence-set ambiguity, finite-sample guarantees, and the price of robustness. No sources found.
7. **Asymptotics in K.** With pairwise information and K → ∞, does the robust value exceed max_i p_i, and how does it scale? With marginals only, Arieli et al. give random dictator. Also open: the correlation gap of log-odds WMRs.
