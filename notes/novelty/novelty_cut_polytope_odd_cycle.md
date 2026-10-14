# Novelty check: odd-cycle / cut-polytope bounds for majority under pairwise co-error bounds

Date of check: 2026-10-13. Literature only; no code was written or run.

Legend: [FT] full text (or full-text extraction) read; [ABS] only abstract, metadata, search listing or a secondary description; [MEM] from memory, not verified in this session; [NA] not accessed.

## 0. Verdict (short)

- Claim (1) (odd-cycle majority bound) is a special case of a classical covering / Boole-type argument. The specific statement for majority voting was NOT found in any source. The pointwise lemma (any m+1 vertices of C_{2m+1} contain an edge, so 1{N>=m+1} <= sum of cycle-edge indicators) is elementary. It is not the sharpest graph-based bound: for symmetric data the fractional covering over all (m+1)-subsets is tighter (see 3.2).
- Claim (2) (bipartite obstruction, randomized rules via cuts) was NOT found in the literature. It reads as a translation of cut-polytope / local-model ideas and needs a rigorous proof.
- Claim (3) (max-min = odd-cycle bound) is open. The right framework (Hailperin LP, quadratic pseudo-Boolean majorants, cut polytope) is known (Hailperin 1965; Boros-Lee 2021). The odd-cycle-only form is probably too narrow. Two concrete reasons: (a) the cut-polytope facet structure changes once the constraint graph (plus apex) has a K5 minor, which for voters means a K4 minor in the constraint graph, not K5; (b) the symmetric case shows the optimal covering is not always a cycle (3.2).
- Overall: PARTIALLY KNOWN. The framework is known. The odd-cycle majority statement looks new as a stated result (not found), but it is a short corollary of known covering ideas. The novelty lies in (i) applying the Hailperin / cut-polytope framework to heterogeneous pairwise co-error bounds for majority juries, and (ii) a correct threshold and refinement of conjecture (3).

## 1. Sources checked

### Cut polytope / K5 minor
- Barahona, F., Mahjoub, A.R. (1986), "On the cut polytope", Math. Programming 36. PDF: https://www.lamsade.dauphine.fr/~mahjoub/Papers/BM86.pdf [FT]. Corollary 3.10: "A graph G is not contractible to K5 if and only if P_c(G) is defined by: 0 <= x(e) <= 1 for each edge e not in a triangle; x(F) - x(C\F) <= |F| - 1 for each chordless cycle C and F subset of C with |F| odd." Theorem 3.3: the odd-cycle inequality is a facet iff C is chordless. Theorem 3.4: 0 <= x(e) <= 1 is a facet iff e is in no triangle. Max-cut is polynomial for graphs not contractible to K5 (cited from Barahona 1983).
- Chimani, Juhnke-Kubitzke, Nover, Roemer, "Cut Polytopes of Minor-Free Graphs", arXiv:1903.01817 [FT]. States the K5-minor-free version ("defined completely by the cycle- and edge inequalities"), attributed to BM86 Sec. 3. Note: the body uses "K5-minor-free" while Barahona (1983) and BM86 say "not contractible to K5". The two are not obviously equivalent. Check before citing the exact hypothesis.
- Barahona, F. (1983), "The max-cut problem on graphs not contractible to K5", Oper. Res. Lett. [ABS, title from Chimani et al.].
- Pitowsky, I. (1991), "Correlation polytopes: their geometry and complexity", Math. Programming 50(3):395-414 [ABS, IAOR abstract]. Defines correlation polytopes, facet families via symmetry, NP-completeness of membership.
- Deza, M., Laurent, M., Geometry of Cuts and Metrics, Springer (1997) [NA]. Source for the correlation-polytope / CUT(n) isomorphism and facet lists.
- Huber, Maric, "Bernoulli correlations and cut polytopes", arXiv:1706.06182 [ABS]. Correlation vectors of n symmetric Bernoulli variables form R(B_n) = 1 - 2 CUT(n) (affine isomorphism); LP can test realizability of a correlation matrix.
- Boros, E., Lee, J. (2021), "Boole's probability bounding problem, linear programming aggregations, and nonnegative quadratic pseudo-Boolean functions", arXiv:2110.10672 [FT]. Hailperin's sharp LP (atomic variables x_S) for P(union) given P(A_i) and P(A_i and A_j) is linked to nonnegative quadratic pseudo-Boolean functions (QPBFs). These correspond to the cut polytope; extremal rays match its facets. Corollaries give LB(H) = LB(QPB), UB(H) = UB(QPB). The paper does not treat "at least k", majority, or odd-cycle bounds in its body.

### Boole / Bonferroni / Hailperin-type bounds
- Hailperin, T. (1965), "Best possible inequalities for the probability of a logical function of events", Amer. Math. Monthly [ABS, title confirmed via a listing; content NA; the general LP for any logical function of events given pairwise info is the key fact]. This is the most important prior art for "exact worst-case probability of a logical function (here majority-failure) given pairwise bounds".
- Prekopa, A. (1988), "Boole-Bonferroni inequalities and linear programming", Oper. Res. 36(1):145-162 [ABS listing].
- Boros, E., Prekopa, A. (1989), "Closed form two-sided bounds for probabilities that at least r and exactly r out of n events occur", Math. Oper. Res. 14(2):317-342 [ABS, title/listing].
- Boros, E., Prekopa, A. (1990), "Sharp bounds on probabilities using linear programming", Oper. Res. 38(2):227-239 [ABS, abstract page]. Sharp LP bounds for P(at least r) and P(exactly r) via dual feasible bases; the abstract does not say whether pairwise information alone suffices.
- Boros, Scozzari, Tardella, Veneziani (2014), "Polynomially computable bounds for the probability of the union of events", Math. Oper. Res. [NA; listing only].
- Veneziani, P. (2008), "Graph-based upper bounds for the probability of the union of events", Electron. J. Combin. 15, #R28 [FT]. LP dual bases; bipartite and threshold-type graph structures appear. Spanning trees and cherry trees appear only in the literature review. No odd-cycle result. No "at least k" result.
- Kounias, E.G. (1968), "Bounds for the probability of a union, with applications", Ann. Math. Statist. 39(6):2154-2158 [ABS via Kuai-Alajaji-Takahara survey, arXiv:1710.07576 (FT of survey)]. Survey states Kounias's bound is Bonferroni-type, maximized over subsets; not cycle-based.
- Hunter (1976), Hunter-type spanning tree upper bound on union [ABS via survey; original NA]. Survey eq. (47): subtracts the largest pairwise sum over a spanning tree.
- Dawson, D., Sankoff, D. (1967), "An inequality for probabilities", Proc. AMS 18:504-507 [ABS via survey]. Survey says it is optimal for bounds using the sum of P(A_i) and the sum of P(A_i and A_j).
- Survey: Kuai, Alajaji, Takahara (2000), IEEE Trans. Inf. Theory 46(7) and "A Short Survey on Bounding the Union Probability using Partial Information", arXiv:1710.07576 [FT]. Covers Kounias, Hunter, Dawson-Sankoff, Galambos-Simonelli (book, 1996). Contains NO odd-cycle bound, NO "at least k" result, NO cut-polytope discussion.
- Chen, J.T. (2003), "A lower bound using Hamilton-type circuits and its applications", J. Appl. Prob. 40(4):1121-1132 [ABS]. Degree-two lower bound on union, with pairwise probabilities as edge weights on a Hamilton-type circuit. A cycle-based bound, but for union lower bounds, not for majority.
- Dohmen, K. (2013), "Lower bounds for the probability of a union via chordal graphs", Electron. Commun. Probab. 18, paper 2357, arXiv:1004.3416 [FT]. Lower Bonferroni-type bound on the union using cliques of chordal graphs, divided by the independence number alpha(G). Says Hunter-type bound is (sum P(A_v) - sum P(A_v and A_w))/alpha(G) for trees. Closest to our "independence number of a graph controls the bound" idea. No odd cycles. Chordality is needed (Remark 2.5).
- Dohmen, K. (2003), Improved Bonferroni Inequalities via Abstract Tubes, Lecture Notes in Math. [ABS, Springer record; volume number MEM].
- Bukszar, J., Prekopa, A. (2001), "Probability bounds with cherry trees", Math. Oper. Res. 26(1):174-189 [ABS via URL and Veneziani's review].
- Kwerel-type "most stringent bounds" for union/intersection of m events under S_1, S_2 (J. Appl. Prob.) [ABS, title confirmed from listing; authors not verified].
- Sathe, Pradhan, Shah (1980), "Inequalities for the probability of the occurrence of at least m out of n events", J. Appl. Prob. [ABS, title via listing; content NA].
- Frolov, A.N. (2020), "On upper and lower bounds for probabilities of combinations of events", arXiv:2008.03588 [FT summary]. Bounds on P(at least r of n) from binomial moments S_1, S_2 (and S_3). Pairwise information only in the simplest case. No graphs, no odd cycles, no LP.
- Adler, Karp, Ross (2025), "Improved Bounds on the Probability of a Union and on the Number of Events that Occur", arXiv:2505.12243 [FT summary]. Bounds on P(X >= r) from S_1..S_{k+1}. Algebraic, no graphs, no odd cycles, no LP.
- Lee, J. (2020), "On the probability of Boolean functions of events in the n-dimensional Euclidean space", arXiv:2001.01368 [FT summary]. Graph of nonempty pairwise intersections is used to reduce inclusion-exclusion. "At least k" covered by LP with S_1..S_3. No cut polytope, no odd cycles.
- Yang, Alajaji, Takahara, "Lower Bounds on the Probability of a Finite Union of Events", arXiv:1401.5543 [FT summary]. LP plus closed forms. No graphs, no at-least-k.
- "Optimal Union Probability Interval Is NP-Hard", arXiv:2605.03556 [NA; title only, text not extractable]. If its title is accurate, exact optimal union-bound computation is NP-hard in general. That fits the view that a closed-form odd-cycle formula could only hold for special structure.

### Majority / Condorcet with dependence (closest in adversary model)
- Kaniovski, S., "Aggregation of correlated votes and Condorcet's jury theorem" (working paper 2008; published Theory and Decision 69(3):453-468, 2010 per listing) [FT of working paper; publication details ABS]. Homogeneous juries: pairwise correlation c is fixed. Min and max of majority correctness over ALL joint distributions consistent with marginals and pairwise correlation are computed by LP. Only n=3 numerics. No graphs, no odd cycles, no cut polytope. This is the closest prior art for the adversary model, but it is homogeneous.
- Zaigraev, A., Kaniovski, S. (2009), "Bounds on the competence of a homogeneous jury", working paper [FT]. Homogeneous exchangeable LP over all higher-order coefficients. Min of majority competence = (n(2p-1)+1)/(n+1), i.e. max P(fail) = 2nq/(n+1) with q = 1-p.
- Kaniovski, S., Zaigraev, A. (2009), "Optimal jury design for homogeneous juries with correlated votes", working paper [FT summary]. Truncated Bahadur representation; not worst case.
- Pivato, M. (2016), "Epistemic democracy with correlated voters", MPRA 69546 [FT summary]. Worst-case average covariance, asymptotic CJT via Chebyshev. No exact worst-case majority bound, no odd cycles.
- Karge, Burkhardt, Rudolph, Rusovac (AAMAS 2024), "To Lead or to be Led: A Generalized Condorcet Jury Theorem under Dependence" [FT summary]. Opinion-leader model, not pairwise worst case.
- Boland, Proschan, Tong (1989), "Modelling dependence in simple and indirect majority systems", J. Appl. Prob. [ABS via citation in Pivato; NA].
- Ladha, K. (1992), "The Condorcet jury theorem, free speech, and correlated votes", AJPS [ABS via citation; NA].
- Majority Correctness in Social Networks, arXiv:2607.14288 [NA; no extractable text].

### Cycle-scenario / Bell-type odd-cycle inequalities (physics literature)
- Araujo, Quintino, Budroni, Cunha, Cabello (2013), "All noncontextuality inequalities for the n-cycle scenario", Phys. Rev. A 88, 022118; arXiv:1206.3212 [ABS; title from listing, abstract fetched concerns complete characterization of classical correlations for consecutive jointly measurable pairs in an n-cycle, generalizing CHSH and KCBS; facet inequalities not verified]. This is the physics analogue of the odd-cycle structure in our claim (1).
- Fine, A. (1982), PRL, "Hidden variables, joint probability, and the Bell inequalities" [MEM].

## 2. Answers to the questions

Q1. Are (1), (2) or (3) stated anywhere?
- Not found in any source. Nothing reviewed states "P(majority of 2m+1 events occur) <= sum of pairwise upper bounds over an odd cycle", or a max-min result for voting with pairwise co-error upper bounds.
- Closest: Kounias/Hunter/Dawson-Sankoff (union), Dohmen 2013 (alpha(G) normalization), Chen 2003 (cycle weights, union lower bound), Veneziani 2008 (graph LP bases), Frolov 2020 and Adler-Karp-Ross 2025 (at least r via binomial moments), Kaniovski 2008 and Zaigraev-Kaniovski 2009 (exact worst-case majority over all joint laws, homogeneous only).

Q2. Known theorem: max P(majority) given P(A_i), P(A_i and A_j) bounds equals an odd-cycle / matching / covering expression?
- Not found. The known exact answer is the Hailperin LP (atomic variables), which equals a QPBF / cut-polytope LP (Boros-Lee 2021). For specific graph-structured answers the known results are Hunter (spanning trees), Veneziani (bipartite and threshold LP bases), Dohmen (chordal cliques with alpha(G)), Chen (Hamilton circuits, lower bound). None addresses at-least-k with odd cycles.
- Fractional covering over (m+1)-subsets is a natural graph-based certificate (3.1). I did not find it stated in the literature for "at least t of n" with pairwise info. The elementary statement 1{|S|>=t} <= C(|S|,2)/C(t,2) is classical in spirit (Boole/Bonferroni family) but I did not verify a source [MEM].

Q3. Cut-polytope: projection onto a graph's edges equals cycle polytope iff no K5 minor?
- Yes for the cut polytope of a graph H (Barahona-Mahjoub 1986 Cor 3.10; Chimani et al. Sec. 2). Statement: P_c(H) = {0<=x(e)<=1 for e not in a triangle; odd chordless cycle inequalities} iff H is not contractible to K5 (BM86), stated as K5-minor-free in Chimani et al. Caveat: contraction vs minor.
- Important correction for our setting. Pairwise moments of K binary errors are the projection of CUT(K+1) (apex 0 with x_{0i} = E_i, and x_{ij} = E_i xor E_j, so E_i E_j = (x_{0i} + x_{0j} - x_{ij})/2). If the constraint graph on voters is G, the feasible pairwise region is P_c(H) with H = G + apex. Then:
  - H has a K5 minor iff G has a K4 minor (apex branch set plus K4 gives K5; a K5 minor in G gives K4 in G). So the odd-cycle description applies iff G is K4-minor-free (series-parallel), not K5-minor-free. With all pairs constrained, K=4 already gives H = K5, with non-cycle facets (pentagonal type of CUT(5); [MEM]: correlation polytope for 4 binary variables has non-triangle facets, the Bell-type pentagonal inequalities).
  - The contraction-vs-minor caveat must be checked for this apex argument.
- Consequence for (3): if G has a K4 minor, extra facets (certificates) exist and may lower the adversary's value below any odd-cycle bound. Whether they actually do is not verified.
- Note: the objective (majority failure) depends on the full law, not only on pairwise moments. The adversary's value equals the concave envelope of the threshold function over the cut vertices of H, maximized over the constraint polytope. Its dual is a min over affine majorants on cut vertices, equivalently QPBF majorants (linear + pairwise + constant). Cycle covers are valid majorants but not all of them.

## 3. My own mathematical notes (not from the literature; to be checked)

### 3.1 Covering certificates
- Pointwise: for a set S of voters with E indicator, 1{|S ∩ E| >= t} <= sum over pairs e in G' with e subset of E of w_e, provided every t-subset T of the voter set has sum_{e subset of T} w_e >= 1 (fractional covering). Then P(fail) <= sum_e w_e * ebar_e.
- Odd cycle C_{2m+1}, t = m+1, weights 1 on cycle edges: this is claim (1). Feasible because any (m+1)-subset of the cycle contains an edge.
- General covering LP: minimize sum_e w_e ebar_e subject to sum_{e in T} w_e >= 1 for all (m+1)-subsets T (of the voter set used by the rule). This is at least as good as the odd-cycle bound for the same rule.

### 3.2 Symmetric (exchangeable) case, closed form
- Symmetric data: P(E_i) = q for all i, P(E_i E_j) <= e for all pairs, n = 2m+1 voters, majority rule (fails if N >= m+1).
- By symmetrization (constraints and objective are permutation-invariant), the adversary can restrict to exchangeable laws, so the exchangeable LP is exact.
- Upper bounds: P(N >= m+1) <= E[N]/(m+1) = nq/(m+1); and P(N >= m+1) <= E[C(N,2)]/C(m+1,2) <= C(n,2) e / C(m+1,2).
- Matching adversary laws (hand-derived): (a) mass a = nq/(m+1) on uniform (m+1)-subsets, rest 0; (b) mass a = C(n,2)e/C(m+1,2) on uniform (m+1)-subsets, rest mass b = nq - a(m+1) on N=1. Each gives symmetric pair co-error, correct marginal q, and the matching value.
- Result: V_sym(n, q, e) = min{ nq/(m+1), C(n,2) e / C(m+1,2) } = (2n/(n+1)) * min{q, 2e}.
- Checks: n = 3: min{3q/2, 3e}. Triangle majority guarantee is 3e when e <= q/2. n = 5: min{5q/3, 10e/3}. The 5-cycle bound of claim (1) gives 5e, which is NOT tight; the complete-graph covering gives 10e/3.
- Consistency: Zaigraev-Kaniovski's homogeneous minimum (n(2p-1)+1)/(n+1) gives max P(fail) = 2nq/(n+1), which equals V_sym in the regime e >= q/2. Good sign for the framework, not a new result.
- Designer side (symmetric case): among majority rules on odd subsets, the smallest (triangle, or a single voter) is best. So the odd-cycle quantity c* = 3e is the relevant one here, and conjecture (3) may hold in the symmetric case restricted to majority rules. It is still unproven for general rules.

### 3.3 Refined conjecture to test instead of (3)
- (3'): For a designer picking an odd subset S and majority on S, the guarantee is min over odd S of the covering LP value of (m+1)-subsets of S with pair weights. Claim (3) would follow if the optimal covering is always a cycle. The symmetric example shows covering with chords can be strictly better for 5 voters.
- (3''): the exact adversary value is the Hailperin / QPBF LP over 2^K atoms; test it numerically on small K (K = 3..6) for heterogeneous ebar and compare with min over odd cycles, min over covering LPs, and (for K >= 4 with complete constraint graph) the pentagonal-facet bounds. This is the decisive test; it needs code, which was out of scope here.

### 3.4 Claim (2) (bipartite)
- The statement "adversary makes one whole side err together" is plausible: it zeroes cross-pair co-errors. Feasibility needs the side marginals to fit (sum of side error probabilities <= 1, and within-side pairs unconstrained). The randomized-rule condition (a distribution over cuts separating each pair with probability >= 1 - ebar/q_min) reads as a local-model / cut-decomposition statement, which is the cut-polytope picture. I did not find a source that states it for voting.

## 4. Most important references to cite

1. Hailperin (1965), "Best possible inequalities for the probability of a logical function of events" [ABS, title] — exact LP for any logical function of events given pairwise information.
2. Boros, Lee (2021), arXiv:2110.10672 [FT] — Hailperin LP and nonnegative QPBF / cut polytope.
3. Barahona, Mahjoub (1986), "On the cut polytope", Math. Programming 36 [FT via PDF] — Cor 3.10, odd-cycle description iff not contractible to K5.
4. Chimani et al., "Cut Polytopes of Minor-Free Graphs", arXiv:1903.01817 [FT] — K5-minor-free version.
5. Barahona (1983), "The max-cut problem on graphs not contractible to K5", Oper. Res. Lett. [ABS, title].
6. Deza, Laurent (1997), Geometry of Cuts and Metrics [NA].
7. Pitowsky (1991), "Correlation polytopes: their geometry and complexity", Math. Programming 50:395-414 [ABS].
8. Huber, Maric, arXiv:1706.06182 [ABS] — R(B_n) = 1 - 2 CUT(n).
9. Prekopa (1988), OR 36:145-162; Boros, Prekopa (1989), MOR 14:317-342; Boros, Prekopa (1990), OR 38:227-239 [ABS] — LP for at least r / sharp bounds.
10. Dohmen (2013), ECP 18, paper 2357 [FT] — chordal, independence-number bounds.
11. Kounias (1968), Ann. Math. Statist. 39:2154-2158; Hunter (1976); Dawson, Sankoff (1967), Proc. AMS 18:504-507 [ABS via survey] — union bounds.
12. Chen (2003), J. Appl. Prob. 40:1121-1132 [ABS] — Hamilton circuits, union lower bound.
13. Veneziani (2008), Electron. J. Combin. 15, #R28 [FT].
14. Frolov (2020), arXiv:2008.03588 [FT summary]; Adler, Karp, Ross (2025), arXiv:2505.12243 [FT summary] — at least r bounds via binomial moments.
15. Zaigraev, Kaniovski (2009), "Bounds on the competence of a homogeneous jury" [FT] and Kaniovski (2008/2010), "Aggregation of correlated votes and Condorcet's jury theorem", Theory and Decision 69:453-468 [FT of working paper] — worst-case majority by LP, homogeneous. These are the closest prior art and must be cited.
16. Pivato (2016), MPRA 69546 [FT summary]; Karge et al. (AAMAS 2024) [FT summary]; Boland, Proschan, Tong (1989) [ABS via citation]; Ladha (1992) [ABS via citation] — jury dependence models.
17. Araujo et al. (2013), PRA 88:022118, arXiv:1206.3212 [ABS] — odd-cycle / n-cycle scenario in physics.

## 5. Could not access (honest list)
- Galambos, Simonelli (1996) book; Kounias (1968) and Hunter (1976) originals (only via survey); Boros-Prekopa 1989/1990 full text (abstracts/listings only); Boros-Scozzari-Tardella-Veneziani (2014) MOR (listing only); Dohmen 2003 LNM (record only); Deza-Laurent book; Pitowsky 1991 full text (abstract only); Huber-Maric (abstract only); Araujo et al. facet statements (abstract only, different content than expected); Sathe-Pradhan-Shah 1980 (listing only); Hailperin 1965 (title only); arXiv:2607.14288 and arXiv:2605.03556 (no extractable text); Ladha 1992 and Boland-Proschan-Tong 1989 (citations only).
- Claims tagged [MEM] (K4/K5 pentagonal facets for 4 binary variables; Bonferroni-type bound P(N>=t) <= C(N,2)/C(t,2) in the literature) were not verified against a source in this session.
