# Novelty check: Boole-type LP, quadratic pseudo-Boolean majorants, self-dual rules, max-cut

Date of check: 2026-10-13. Scope: probability bounds for Boolean functions of binary errors X in {0,1}^K with known P(X_i=1)=q_i and pairwise constraints; LP duality via quadratic pseudo-Boolean majorants; claims A4, A5, A8.

Tags: [FT] = body-level text obtained (via fetch extraction of the full PDF/HTML, not a line-by-line proof check); [ABS] = abstract or metadata page only; [NA] = not accessible, or only a bibliographic record / secondary citation was found. Year and volume entries marked "inferred" come from DOI or volume numbering and were not confirmed on a page.

---

## 1. Source list

### Core: Boole / Hailperin LP and quadratic majorants

1. **Hailperin, T. (1965).** "Best possible inequalities for the probability of a logical function of events." *American Mathematical Monthly* 72, pp. 343-359. DOI 10.2307/2313491. Tag: [NA] (metadata only; MaRDI and zbMATH records confirmed the title, venue, year and DOI; content not read). Cited as the origin of the LP formulation by Boros-Lee and Kaski et al.

2. **Hailperin, T. (1996).** *Sentential Probability Logic: Origins, Development, Current Status, and Technical Applications.* Lehigh University Press. Tag: [ABS] (publisher page confirms the title and year; content not read). Also *Boole's Logic and Probability* (North-Holland, 1986). [ABS] (publisher listing only).

3. **Boros, E. and Lee, J. (2026).** "Boole's probability bounding problem, linear programming aggregations, and nonnegative quadratic pseudo-Boolean functions." *Mathematics of Operations Research* 51(1), pp. 134-148. DOI 10.1287/moor.2023.0019. arXiv:2110.10672 (v4, Jan 2025). Tag: [FT] (arXiv HTML extraction: sections, Theorems 1-8, Corollaries 1-3, and quoted statements).
   - What it states: Hailperin's LP for Boole's problem is hard and has exponentially many variables. Bounds come from aggregations of the LP. The image of the nonnegative orthant under the aggregation is a cone of QPBF-indexed inequalities (Theorem 8). Corollary 2: LB(H)=LB(QPB) and UB(H)=UB(QPB), i.e., the Hailperin optimum equals the optimum over the nonnegative-QPBF cone. Corollary 3 compares several aggregation bounds.
   - **Critical differences from our setting (verified in the text):**
     (a) Boole's problem is posed with **equalities** P(A_i)=p_i and P(A_i∩A_j)=p_ij. Our setting uses **inequalities** P(X_i=X_j=1) <= b_ij. Under inequalities the pair multipliers gamma_ij must be >=0. Under equalities they are free. So our dual is not the Boros-Lee dual as stated.
     (b) "Nonnegative QPBF" in Boros-Lee means nonnegative **on the cube** (F_alpha(Z)>=0 for all Z in {0,1}^V), not nonnegative coefficients.
     (c) The paper does **not** present the bound as minimizing a + sum q_i beta_i + sum b_ij gamma_ij over majorants Q >= F. It uses a primal LP over aggregated variables y.
     (d) The paper does **not** discuss self-dual functions, threshold functions, inequality constraints on pairs, or a condition sum q <= 1.
     (e) Cut-polytope equivalence: the paper says the cone of nonnegative QPBFs "is known to be equivalent with the cut-polytope" and cites **C. De Simone, "The cut polytope and the boolean quadric polytope"** (not verified; not accessed). It does not cite Boros-Hammer 1993 for this.

4. **Boros, E. and Prékopa, A. (1989).** "Closed form two-sided bounds for probabilities that at least r and exactly r out of n events occur." *Mathematics of Operations Research* 14(2), pp. 317-342 (DOI 10.1287/moor.14.2.317; year inferred from volume). Tag: [ABS] (IAOR and INFORMS abstracts only).
   - States: closed-form lower and upper bounds for P(at least r) and P(exactly r) using binomial moments; monotone decrease of dual vector components for P(at least one).
   - Implies A5 (the bound P(N>=3) <= 1 - S1 + S2)? **Partially, unverified.** The abstract does not state the r=3 case or the S1,S2-only form. The bound is a valid special case of Bonferroni-type inequalities (see Section 2). Full text not read.

5. **Prékopa, A. (1988).** "Boole-Bonferroni inequalities and linear programming." *Operations Research* 36(1), pp. 145-162. DOI 10.1287/opre.36.1.145. Tag: [ABS].
   - States: sharp lower and upper bounds on P(at least one event) via LP over known binomial moments, using only the first few terms of inclusion-exclusion.
   - Implies A4/A5/A8? **No.** Covers only the union (at least one) case with binomial-moment inputs. Does not treat general Boolean functions or self-dual rules.

6. **Prékopa, A. (1990).** "Sharp bounds on probabilities using linear programming." *Operations Research* 38(2), pp. 227-239. Tag: [ABS].
   - States: sharp bounds via LP over binomial moments for **exactly r** and **at least r** events, based on dual feasible basis characterization theorems; also linear inequalities among binomial moments.
   - Implies A5? **Partially.** The LP framework over S_j for at-least-r is the closest prior art for the A5 counterexample. The specific use of inequality constraints b_ij (not equalities) and the conclusion about sum q > 1 are not addressed in the abstract. Does not address self-dual or graph-based min-max.

7. **Kwerel, S. M. (1975).** "Most stringent bounds on the probability of the union and intersection of m events for systems partially specified by S_1, S_2, ..., S_k, 2<=k<m." *Journal of Applied Probability* 12(3), pp. 612-619. DOI 10.2307/3212879. Tag: [ABS].
   - States: the most stringent upper and lower bounds on P(union) and P(all events) given S_1,...,S_k (binomial moments).
   - Implies A5? **Partially.** Same family as above (union and intersection, not at-least-3). Confirms that "most stringent bounds given S_1,S_2" is classical.

8. **Padberg, M. (1989).** "The Boolean quadric polytope: some characteristics, facets and relatives." *Mathematical Programming* 45(1), pp. 139-172. Tag: [ABS].
   - States: polyhedral study of the Boolean quadric polytope (linearized quadratic 0/1 form); facet families; maximizing nonnegative quadratic forms by LP rounding.
   - Implies A8? **Partially.** Relevant to the LP side of the quadratic majorant dual (BQP as the convex hull of (x, xx^T)). The cut polytope is not mentioned in the abstract.

9. **Boros, E. and Hammer, P. L. (1993).** "Cut-polytopes, boolean quadric polytopes and nonnegative quadratic pseudo-Boolean functions." *Mathematics of Operations Research* 18(1), pp. 245-253 (DOI 10.1287/moor.18.1.245; year inferred from volume). Tag: [ABS].
   - States: one-to-one correspondence between valid inequalities (facets) of the cut polytope of K_{n+1} and nonnegative extremal quadratic pseudo-Boolean functions in n variables; large facet classes; polynomial separation.
   - Implies A8? **Partially.** Establishes the cut-polytope <-> nonnegative QPBF correspondence (the same cone Boros-Lee use). It does not address a min over graph bipartitions, a max-cut separation for a Boole-type LP, or pair-inequality constraints.

10. **Pitowsky, I. (1991).** "Correlation polytopes: their geometry and complexity." *Mathematical Programming* 50(3), pp. 395-414. Tag: [ABS].
    - States: correlation polytopes arising in probability and propositional logic; symmetry and facet families; membership is NP-complete.
    - Implies A8? **Partially.** Gives the complexity context for the correlation (cut-type) polytope side. It does not mention cut polytope equivalence or Boole's problem on the page.

11. **Boros, E. and Hammer, P. L. (2002).** "Pseudo-Boolean optimization." *Discrete Applied Mathematics* 123(1-3), pp. 155-225. DOI 10.1016/S0166-218X(01)00341-9. Tag: [NA] (metadata only; content and roof-duality / posiform details not verified).

12. **Crama, Y. and Hammer, P. L. (2011).** *Boolean Functions: Theory, Algorithms, and Applications.* Cambridge University Press. Tag: [ABS] (table of contents only). Chapter 1 is "Fundamental concepts and applications". Part II special classes (not listed by title) includes threshold functions per the Chapter 1 summary. Chapter 4 is "Duality theory". Self-dual functions are not named in the visible text. Posiforms are not mentioned. **Self-dual and threshold characterizations were not verified from this source.**

### Union bounds, graph and tree structures

13. **Veneziani, P. (2008).** "Graph-based upper bounds for the probability of the union of events." *The Electronic Journal of Combinatorics* 15, #R28. Tag: [FT] (PDF extraction; propositions and citations).
    - States: LP-dual feasible-basis bounds on P(union) from individual and pairwise probabilities. Built on bipartite partitions (A,B): Prop. 3 and Prop. 5 give bounds with negative cross-part terms sum_{i in A, j in B} p_ij. Prop. 7 is a nested "threshold-type" structure. Cites Hunter, Prékopa et al., Bukszár-Prékopa, Bukszár-Szántai, Kwerel, Dawson-Sankoff, Hailperin.
    - Implies A4 or A8? **Partially.** The cross-part (bipartition) term resembles the cut term in B_graph. Each bound is for a fixed partition; no min over partitions of a max of within-part sums; no self-dual or threshold Boolean functions; not framed as max-cut. Claims dual-basis sharpness only for its numerical example.

14. **Boros, E., Scozzari, A., Tardella, F., Veneziani, P. (2014/2015).** "Polynomially computable bounds for the probability of the union of events." *Mathematics of Operations Research*. DOI 10.1287/moor.2014.0657. Volume and year not shown on the page; not verified. Tag: [ABS] (INFORMS abstract, paraphrased).
    - States: best bounds come from exponential-size LPs; the authors add constraints to dual LPs to get three polynomially solvable bound classes; two bounds coincide with classical bounds at polyhedron vertices.
    - Implies A4/A8? **No** for self-dual or graph min-max. **Partially** in the sense that dual-LP tightening for Boole's union bound is known.

15. **Bukszár, J. and Prékopa, A. (2001).** "Probability bounds with cherry trees." *Mathematics of Operations Research* 26(1), pp. 174-187 (start page from DOI moor.26.1.174; end page and year-of-volume inferred). Tag: [ABS].
    - States: third-order upper bounds on P(union) via cherry trees; any Hunter-Worsley bound is majorized by a t-cherry-tree bound; cherry-tree bounds are feasible solutions of the dual of the Boolean probability bounding problem.
    - Implies A4 or A8? **Partially.** Gives a tree-graph dual-LP representation for union bounds from pairwise information. Does not treat self-dual rules or max-cut.

16. **Hunter, D. (1976)** upper bound via spanning trees; **Kounias (1968)** and **Dawson, D. A. and Sankoff, D. (1967)** bounds on the union. Tag: [NA] (not accessed directly; cited via Veneziani 2008 and Bukszár-Prékopa).

17. **Bukszár, J., Mádi-Nagy, G., Szántai, T. (2012).** "Computing bounds for the probability of the union of events by different methods." *Annals of Operations Research* 201, pp. 63-81. DOI 10.1007/s10479-012-1231-1. Tag: [ABS]. Compares hypermultitree and binomial-moment methods numerically. Implies A4/A8? **No.**

18. **Yang, J., Alajaji, F., Takahara, G. (2017).** "A short survey on bounding the union probability using partial information." arXiv:1710.07576. Tag: [ABS] (abstract only). Covers pairwise-information bounds and the Gallot-Kounias bound. Does not mention Boros-Lee, QPBFs, max-cut or self-dual functions. Implies A4/A5/A8? **No** (no sharp-LP or graph-min-max result on the abstract).

19. **Padmanabhan, D. and Natarajan, K. (2019).** "Tree bounds for sums of Bernoulli random variables: a linear optimization approach." arXiv:1910.06321. Tag: [ABS]. Tightest bounds on P(sum > k) given pairwise distributions. NP-hard on complete graphs; polynomial on trees via LP. Implies A8? **Partially** (complexity of pairwise-specified bounds; the complete-graph hardness is the analogue of the Kaski et al. result). Does not treat QPBFs, max-cut or self-dual functions.

20. **Ramachandra, A. and Natarajan, K. (2020/2022).** "Tight probability bounds with pairwise independence." arXiv:2006.00516 (v6). Tag: [ABS]. Closed-form tightest union bound under pairwise independence (ratio to Boole bound at most 4/3). No LP duality, self-dual or max-cut content on the page. Implies A4/A5/A8? **No**, though relevant as "pairwise-determined tight bound" precedent.

21. **Frolov, A. N. (2020).** "On upper and lower bounds for probabilities of combinations of events." arXiv:2008.03588. Tag: [ABS]. Bounds on P(exactly r) and P(at least r) that can be equalities. No LP, self-dual or max-cut content. Implies A5? **Partially** (at-least-r bounds, but not from pairwise data and not Fréchet-optimal in the abstract).

### Hardness and cut structure

22. **Kaski, P., Mannila, H., Mohapatra, C. K. (2026).** "Optimal union probability interval is NP-hard." ESA 2026 (34th European Symposium on Algorithms), LIPIcs vol. 388, article 129, pp. 129:1-129:19. DOI 10.4230/LIPIcs.ESA.2026.129. arXiv:2605.03556. Tag: [FT] (PDF extraction; Theorems 2, 5, 6 summarized).
    - States: Hailperin's LP geometry via coordinate projections of an elementary polytope on Venn atoms; computing the min or max union probability given feasible **pairwise** intersection probabilities is NP-hard (Thm 2), even with all sets of size <= 2 (reduction from fractional graph coloring for min, k-clique for max). Membership in correlation polytope NP-hard (Thm 5); union polytope NP-hard (Thm 6). Cites Hailperin, Boros-Lee, cut and correlation polytopes.
    - Implies A8? **Partially.** Gives hardness of the exact union interval from pairwise information, which is the setting in which B_graph is an approximation. It does **not** treat self-dual functions, threshold functions, inequality pair constraints, or min over Boolean rules. Also note: its min/max is over probabilities with equalities (as in Boole's problem), not with upper bounds only, as far as the extraction shows.

23. **Grötschel, M. and Pulleyblank, W. R. (1981).** "Weakly bipartite graphs and the max-cut problem." *Operations Research Letters* 1(1), pp. 23-27. Tag: [NA] (BibTeX metadata only; content not read). Relevance: LP relaxations of max-cut and odd cycles, the natural comparator for our LP whose separation is max-cut.

24. **Kral, Sereni, Stacho (2011).** "Min-max relations for odd cycles in planar graphs." arXiv:1108.4281. Tag: [ABS]. Vertex-based odd-cycle packing and transversal; not edge-weighted, not max-cut LP. Implies A8? **No** (only a peripheral reference for the odd-cycle transversal min-max).

### Min-max over rules and voting under dependence

25. **Arieli, I., Babichenko, Y., Talgam-Cohen, I., Zabarnyi, K. (2023).** "A random dictator is all you need." arXiv:2302.03667 (v2 Nov 2023). Tag: [ABS].
    - States: symmetric agents with known marginal recommendation distributions, adversarial correlation among recommendations; under minimax, regret and approximation-ratio criteria the optimal rule for many agents is the random dictator.
    - Implies A4 or A5? **Partially, important caveat.** It shows that the minimax optimum over aggregation rules under adversarial correlation can be a **randomized** rule. This matters for our A4 (min over deterministic self-dual F). See Section 3 point 3.

26. **Arieli, I., Babichenko, Y., Smorodinsky, R. (2017/2018).** "Robust forecast aggregation." arXiv:1710.02838. Tag: [ABS]. Regret-based aggregation of expert forecasts for binary states. No Boolean-function, self-dual or LP content on the page. Implies A4/A5/A8? **No.**

27. **Frongillo, R., Monroe, M., Neyman, E., Waggoner, B. (2025).** "Robust forecast aggregation via additional queries." arXiv:2512.05271. Tag: [ABS]. Not about Boolean rules under pairwise constraints. Implies? **No.**

28. **Kaniovski, S. and Zaigraev, A. (2009).** "Optimal jury design for homogeneous juries with correlated votes." Tag: [FT] (extraction). Homogeneous jury, common pairwise correlation; exact correctness under specified correlation (Bahadur representation); quota rules. No worst-case (Fréchet/Boole) bounds, no LP, and no self-dual analysis. Implies A4/A5? **No.** Note the paper observes that P(X_i=X_j=1) <= min(q_i,q_j), a Fréchet constraint used only for feasibility.

29. **Karge, J., Burkhardt, J.-M., Rudolph, S., Rusovac, D. (2024).** "To lead or to be led: a generalized Condorcet jury theorem under dependence." AAMAS 2024. Tag: [ABS]. Opinion-leader dependence; no pairwise worst-case analysis. Implies? **No.**

30. **Ruta, D. and Gabrys, B. (2002).** "A theoretical analysis of the limits of majority voting errors for multiple classifier systems." *Pattern Analysis and Applications* 5(4), pp. 333-350. Tag: [ABS]. Discrete error distributions; no pairwise Boole-type bound stated. Implies? **No.**

31. **Kuncheva, L. I., Whitaker, C. J., Shipp, C. A., Duin, R. P. W. "Is independence good for combining classifiers?"** (c. 2000). Tag: [FT] (extraction). Empirical and enumerated min/max majority-vote accuracy with pairwise kappa; no analytical bounds; no LP; no sharpness proof. Implies? **No.**

32. **Zhu, M. (2018).** "When is the majority-vote classifier beneficial?" arXiv:1307.6522. Tag: [ABS]. CLT analysis under pairwise error correlation; no Boole/Fréchet bounds, no self-dual functions. Implies? **No.**

33. **Rohwer, R. (2020).** "The optimal 'AND'." arXiv:2005.11635. Tag: [ABS]. Bayesian (Jeffreys prior, KL loss), not worst-case. Implies? **No.**

34. **Rüger, B. (1978).** "Das maximale Signifikanzniveau des Tests ... " *Metrika* 25, pp. 171-178. DOI 10.1007/BF02204362. Tag: [ABS]. Smallest upper bound on P(at least k of n events) given individual probabilities only (first-order), and equal to the bound for "exactly k". Implies A5? **No** (marginal-only, no pairwise terms). Useful as a classical min-max over "at least k" with only marginal information.

35. **Pitowsky, I. (1989).** *Quantum Probability - Quantum Logic.* Lecture Notes in Physics 321 (Springer). Tag: [NA] (cited by Kaski et al. for the correlation-polytope / Boole problem gap; not accessed).

### Background (metadata only)

36. **Galambos, J. and Simonelli, I. (1996).** *Bonferroni-type Inequalities with Applications.* Springer. Tag: [NA] (publisher page only; content not read). Standard reference for generalized Bonferroni bounds on P(at least r).

37. **Prékopa, A. (publication list, rutcor.rutgers.edu).** Tag: [ABS] (title list only). Confirms the titles and years of the Boole-Bonferroni, Boros-Prékopa 1989, Sharp-LP 1990, Discrete Moment Problems, cherry-trees 2001 and union-bound papers.

---

## 2. Claim-by-claim check

### A5: P(N>=3) <= 1 - S1 + S2 (and the failure of the formula for sum q > 1)

Our own verification:
- For integer N, C(N-1,2) >= 1{N>=3}: values are 1 (N=0), 0 (N=1), 0 (N=2), 1 (N=3), 3 (N=4), 6 (N=5). Hence P(N>=3) <= E[C(N-1,2)] = 1 - S1 + S2 with S_j = E[C(N,j)]. This matches the Q in the dual: Q(x)=1 - sum x_i + sum_{i<j} x_i x_j.
- With S1=1.5, S2 <= 0.5 (pair bounds 0.05 times 10 pairs): 1 - S1 + S2 <= 0, so P(N>=3)=0 for every feasible law. The LP value for F = "at least 3 of 5" is 0, which is below B_graph = 0.15.
- Majority of 5 (odd K) is self-dual, so it is in the min over self-dual F. Its min-max value is 0, so the formula min{q_min, B_graph}=0.15 fails for sum q > 1, as claimed.
- Feasibility check: N in {1,2} with P(N=1)=P(N=2)=0.5 gives S1=1.5, S2=0.5, and symmetric pair probabilities 0.05 are realizable. The law exists, so the counterexample is not vacuous.

Literature:
- The bound is a special case of the classical Bonferroni / binomial-moment family (Prékopa 1988, 1990; Boros-Prékopa 1989; Kwerel 1975; Galambos-Simonelli 1996). Whether the specific form 1 - S1 + S2 for at-least-3 appears there could not be verified (full texts not accessed).
- The counterexample itself (the failure of a min-max formula when sum q > 1) is not in any source found.

Verdict A5: The inequality and its technique are known in substance (partially known). The specific counterexample to the formula is our observation; not found in the literature.

### A4: min over self-dual F of the min-max equals min{q_min, B_graph} when sum q <= 1

Our spot checks (consistent with A4):
- K=2: self-dual F are the dictators and their complements. The min-max is q_min. B_graph = infinity (g(y)=0 since no bipartition has an internal edge on both sides), so the formula gives q_min. Consistent.
- K=3 triangle, q=0.2, b=0.1: dictators give 0.2. Majority gives at most P(N>=2) <= S2 = 0.3 and <= S1/2 = 0.3. Here g(y)=min_e y_e and B_graph = 0.3, so the formula gives min{0.2, 0.3} = 0.2. Consistent.

Reformulation (our derivation, checked by hand; may help the proof):
- Using LP duality for each fixed F, min over F of max over laws = min over Q of val(Q), where the minimum runs over Q with Q(x) >= F(x) for some self-dual F.
- For a fixed Q, such an F exists iff Q(x) >= 0 for all x in {0,1}^K and, for each complementary pair {x, 1-x}, max(Q(x), Q(1-x)) >= 1.
- So the self-dual min-max equals min val(Q) over quadratic Q with gamma >= 0, Q >= 0 on the cube, and the covering condition above. This is a cleaner statement than the Boros-Lee dual, and it is not in the literature we found.

Important caveat (not in the user's note):
- The claim is for **deterministic** self-dual F. Min over mixed rules (randomized aggregation) is a bilinear game between rules and laws, so by von Neumann minimax the randomized value equals max over laws of min over self-dual F. This can be strictly smaller than the deterministic min-max. Arieli et al. (2023) show randomization is essential in a related minimax aggregation problem. Please check whether A4 is meant for deterministic F, and state it.
- The normalization step in the proof sketch ("zero constant term, nonnegative linear terms, possible iff sum q <= 1") was not checked against literature and was not verified here.

Literature:
- Duality of Hailperin's LP with quadratic majorants (Hailperin 1965; Boros-Lee 2026, equality form, Corollary 2) is known.
- Graph and bipartition-based pairwise bounds for union (Veneziani 2008, Props 3, 5, 7; Hunter; Bukszár-Prékopa 2001) are known.
- Min over self-dual or threshold rules of max error under pairwise constraints: not found.
- B_graph as defined (inf over y >= 0 of b.y / g(y), with g the min-max within-part weight): not found.

Verdict A4: Not found in the stated form (min over self-dual F, inequality pair bounds, B_graph). The underlying duality is known in equality form (Boros-Lee 2026), and graph/bipartition-based union bounds are known (Veneziani 2008). Treat as **not found**, with partial prior art listed. Add the deterministic-vs-randomized caveat.

### A8: F(y)/2 <= g(y) <= F(y), with F(y) = y(E) - MaxCut(y), so B_graph is within a factor 2 of inf b.y/F(y), an LP with max-cut separation

Our verification:
- Let T = y(E). For a bipartition S, y(E(S)) + y(E(S^c)) = T - cut_y(S), so F(y) = min_S [y(E(S)) + y(E(S^c))] = T - MaxCut(y).
- g(y) = min_S max{a_S, b_S} with a_S + b_S = T - cut_y(S). Since max{a,b} <= a+b, g <= F. Since max{a,b} >= (a+b)/2, g >= F/2. Both hold. Consequently B_graph lies in [inf b.y/F, 2 inf b.y/F].
- The LP min b.y s.t. F(y) >= 1, y >= 0 is equivalent to min b.y s.t. y(E) - cut_y(S) >= 1 for all S, which has max-cut as the separation oracle. Equivalently, y(T') >= 1 for every edge set T' whose removal leaves a bipartite graph (odd-cycle transversal). So the LP is a blocker-type covering LP over odd-cycle transversals. (Our observation; not checked against the weakly-bipartite or blocker literature in detail.)

Literature:
- The cut-polytope <-> nonnegative QPBF correspondence is known (Boros-Hammer 1993; De Simone 1990 as cited by Boros-Lee 2026).
- Correlation-polytope membership is NP-complete (Pitowsky 1991). The exact union interval from pairwise information is NP-hard (Kaski-Mannila-Mohapatra, ESA 2026).
- Max-cut LP relaxations for weakly bipartite graphs: Grötschel-Pulleyblank (1981; metadata only).
- The specific LP min b.y s.t. F(y) >= 1 (with F the bipartization function) as a bound for Boole-type probabilities: not found.

Verdict A8: The inequality F/2 <= g <= F is our elementary verification. Max-cut and cut-polytope links to Boole-type bounds are **partially known** (Boros-Hammer 1993; Boros-Lee 2026 via De Simone; Pitowsky 1991; Kaski et al. 2026). The LP with max-cut separation as a Boole-type bound is **not found**.

---

## 3. Other findings that change the framing

1. **Equality vs inequality.** Boros-Lee (2026) and Kaski et al. (2026) both treat pairwise probabilities as equalities. Our b_ij constraints are inequalities, which changes signs of dual variables (gamma >= 0). A reviewer could say "this is Boole's problem with a different constraint set." Make the distinction explicit in the paper.
2. **Nonnegative QPBF.** Boros-Lee's "nonnegative" means nonnegative on the cube. Our Q satisfies Q >= F >= 0 on the cube, so our majorants are nonnegative QPBFs in their sense, plus gamma >= 0. Use the same terminology and cite Boros-Lee.
3. **Randomized rules.** See the caveat under A4. The minimax over randomized self-dual rules equals max over laws of min over self-dual F. Arieli et al. (2023, random dictator) show randomization can be the minimax optimum under adversarial correlation with marginal-only information.
4. **Hardness.** Kaski et al. (ESA 2026) show that the exact union interval from pairwise information is NP-hard. This supports presenting B_graph as a polynomial-time surrogate with a factor-2 guarantee. Please confirm that the hardness covers the exact interval for our inequality setting, which the extraction does not show.
5. **Prior art for threshold-type optimal rules under pairwise ambiguity.** None found. The closest are Arieli et al. (2023) (marginals only, random dictator), Kaniovski-Zaigraev (2009) (no worst case), and Veneziani (2008) (threshold-type graph structure, but for union bounds).

---

## 4. Verdicts

| Claim | Verdict | Closest known statement |
|---|---|---|
| A4: min over self-dual F of min-max = min{q_min, B_graph} for sum q <= 1 | **Not found** (stated form). Partially known components. | Boros-Lee 2026 Cor. 2 (Hailperin LP = nonnegative-QPBF dual, equality form); Veneziani 2008 (bipartition-based union bounds); Bukszár-Prékopa 2001 (cherry-tree dual feasible bases). Deterministic vs randomized caveat applies. |
| A5: P(N>=3) <= 1 - S1 + S2; counterexample to formula for sum q > 1 | **Partially known.** Bound is classical in substance; counterexample is ours. | Prékopa 1988/1990; Boros-Prékopa 1989; Kwerel 1975; Galambos-Simonelli 1996 (binomial-moment / Bonferroni family). Exact at-least-3 form not confirmed in accessible text. |
| A8: F/2 <= g <= F, factor-2 approximation with max-cut separation | **Inequality: ours (elementary, verified). LP with max-cut separation as Boole-type bound: not found.** Cut polytope link partially known. | Boros-Hammer 1993 (cut polytope <-> nonnegative QPBF); Boros-Lee 2026 (citing De Simone 1990); Pitowsky 1991 (correlation polytope membership NP-complete); Kaski-Mannila-Mohapatra 2026 (NP-hard union interval). |

---

## 5. What we must cite

Core (must cite):
- Hailperin, T. (1965). Amer. Math. Monthly 72, 343-359. DOI 10.2307/2313491. [NA]
- Hailperin, T. (1996). Sentential Probability Logic. Lehigh UP. [ABS]
- Boros, E., Lee, J. (2026). Boole's probability bounding problem, LP aggregations, and nonnegative QPBFs. MOR 51(1), 134-148. arXiv:2110.10672. [FT]
- Boros, E., Hammer, P. L. (1993). Cut-polytopes, boolean quadric polytopes and nonnegative quadratic pseudo-Boolean functions. MOR 18(1), 245-253. [ABS]
- De Simone, C. (1990). The cut polytope and the boolean quadric polytope. (Cited in Boros-Lee; venue and year not verified.) [NA]
- Boros, E., Hammer, P. L. (2002). Pseudo-Boolean optimization. DAM 123(1-3), 155-225. [NA]
- Kaski, P., Mannila, H., Mohapatra, C. K. (2026). Optimal union probability interval is NP-hard. ESA 2026, LIPIcs 388, 129. arXiv:2605.03556. [FT]

Binomial-moment and Boole-Bonferroni LP (must cite for A5):
- Prékopa, A. (1988). Boole-Bonferroni inequalities and linear programming. OR 36(1), 145-162. [ABS]
- Boros, E., Prékopa, A. (1989). Closed form two-sided bounds for probabilities that at least r and exactly r out of n events occur. MOR 14(2), 317-342. [ABS]
- Prékopa, A. (1990). Sharp bounds on probabilities using linear programming. OR 38(2), 227-239. [ABS]
- Kwerel, S. M. (1975). Most stringent bounds on the probability of the union and intersection of m events. JAP 12(3), 612-619. [ABS]
- Galambos, J., Simonelli, I. (1996). Bonferroni-type Inequalities with Applications. Springer. [NA]

Graph and tree bounds (must cite for A4 and A8 context):
- Veneziani, P. (2008). Graph-based upper bounds for the probability of the union of events. EJC 15, #R28. [FT]
- Boros, E., Scozzari, A., Tardella, F., Veneziani, P. (MOR, DOI 10.1287/moor.2014.0657). [ABS]
- Bukszár, J., Prékopa, A. (2001). Probability bounds with cherry trees. MOR 26(1), 174-187. [ABS]
- Hunter, D. (1976); Kounias (1968); Dawson, Sankoff (1967). Cited via Veneziani 2008. [NA]

Cut polytope and hardness:
- Padberg, M. (1989). The Boolean quadric polytope. Math. Prog. 45(1), 139-172. [ABS]
- Pitowsky, I. (1991). Correlation polytopes: their geometry and complexity. Math. Prog. 50(3), 395-414. [ABS]
- Grötschel, M., Pulleyblank, W. R. (1981). Weakly bipartite graphs and the max-cut problem. ORL 1(1), 23-27. [NA]
- Crama, Y., Hammer, P. L. (2011). Boolean Functions. CUP. [ABS]

Robust aggregation and voting (cite as the nearest prior art):
- Arieli, I., Babichenko, Y., Talgam-Cohen, I., Zabarnyi, K. (2023). A random dictator is all you need. arXiv:2302.03667. [ABS]
- Kaniovski, S., Zaigraev, A. (2009). Optimal jury design for homogeneous juries with correlated votes. [FT]
- Ruta, D., Gabrys, B. (2002). Limits of majority voting errors for multiple classifier systems. PAA 5(4), 333-350. [ABS]

---

## 6. Recommended wording for the paper

- Position the LP as an inequality-constrained variant of the Hailperin/Boole LP, citing Boros-Lee (2026) for the equality case and the QPBF duality.
- Present A5 as a Bonferroni-type bound from the binomial-moment family (Prékopa; Kwerel), with the counterexample as our observation.
- Present A8 as an elementary factor-2 comparison, and cite the cut polytope and correlation polytope literature (Boros-Hammer 1993; Pitowsky 1991; Kaski et al. 2026) for the complexity context.
- For A4, state explicitly that the min is over deterministic self-dual rules, and either prove the randomized version or remark on the gap (Arieli et al. 2023).
- Before submission: obtain full texts of Boros-Prékopa 1989, Prékopa 1990 and De Simone 1990 to confirm whether the at-least-3 bound 1 - S1 + S2 and the cut-polytope correspondence appear as stated; these were not accessible.

## 7. Access log

- Full text or detailed extraction obtained: Boros-Lee 2026 (arXiv HTML), Kaski et al. 2026 (PDF), Veneziani 2008 (PDF), Kaniovski-Zaigraev 2009 (PDF), Kuncheva et al. (PDF), Prékopa publication list (site).
- Abstract or metadata only: Boros-Prékopa 1989, Prékopa 1988 and 1990, Kwerel 1975, Boros-Hammer 1993, Boros-Scozzari-Tardella-Veneziani, Padberg 1989, Pitowsky 1991, Crama-Hammer 2011 (TOC), Bukszár-Prékopa 2001, Arieli et al. 2023, Frolov 2020, Ramachandra-Natarajan, Padmanabhan-Natarajan, Yang-Alajaji-Takahara 2017, Ruta-Gabrys 2002, Zhu 2018, Rohwer 2020, Frongillo et al. 2025, Karge et al. 2024, Rüger 1978.
- Not accessible: Hailperin 1965 (content), Hailperin 1996 (book), Boros-Hammer 2002 (content), De Simone 1990, Grötschel-Pulleyblank 1981 (content), Galambos-Simonelli 1996 (content), Hunter 1976, Dawson-Sankoff 1967, Pitowsky 1989 (book), arXiv 1710.07576 (full text not read).
- Fetch errors: first attempt at arXiv abs failed on URL scheme (retried with https). Springer link s10479-014-1712-5 was a different paper (Unuvar-Özguven-Prékopa 2015), not Boros-Scozzari-Tardella-Veneziani; the correct record was found via INFORMS DOI.
