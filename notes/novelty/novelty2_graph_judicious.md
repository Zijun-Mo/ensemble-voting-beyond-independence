# Novelty check 2: judicious bipartitions, bipartization LP, B_graph vs odd cycles

Date of check: 2026-10-13. Literature search only, no code.

Tags: [FT] = PDF fetched and body extracted (extraction is a summarizing tool, not verbatim; re-check constants against originals before quoting). [ABS] = abstract or summary page only. [NA] = not accessible, or cited only from another source, or from memory.

## 0. Verified identities (my derivations, checked by hand)

Notation: y weights on K_K, F(y) = min_S internal(S) = y(E) - MaxCut(y), g(y) = min_S max{internal(S), internal(V\S)}.
- g(y) >= F(y)/2: max(e1,e2) >= (e1+e2)/2 >= F/2.
- g(y) <= F(y): take S = a max-cut partition, then e1+e2 = F so max <= F.
- g(y) <= y(E)/3 for every weighted y: at a max cut, 2e_i <= cut (local optimality), so M = max(e1,e2) <= cut/2 and y(E) = e1+e2+cut >= 3M. Equality for K3. This is the weighted form of the survey's g(m) <= m/3.
- A3 proof: a wrong-vote set W satisfies 2y(E(W)) + y(delta(W)) > y(E), which gives y(E(W)) > y(E(V\W)), so y(E(W)) >= g. Markov over E[y(E(W))] = sum y_ij P(i,j both wrong) <= sum b_ij y_ij gives A3. Verified.
- A8 sandwich: g <= F gives B_bip <= B_graph; g >= F/2 gives B_graph <= 2 B_bip. Verified.
- B_bip as an LP: constraint "internal(S) >= 1 for all S" is equivalent to F(y) >= 1, and by homogeneity B_bip = inf b.y / F(y). Separation problem = min_S internal(S) = F(y) = weighted max-cut (NP-hard in general).
- LP dual (standard LP duality, my derivation): B_bip = max over probability distributions mu on bipartitions of min_e b_e / P_mu(e is internal).
- Integer version: a 0/1 y is feasible iff its support is non-bipartite, so integer B_bip = shortest odd cycle (BFS-computable).
- Uniform b (b_ij = 1): B_graph = 3 (triangle, using g <= y(E)/3). B_bip = 2 (inf, approached by large complete graphs since maxcut(K_n)/m -> 1/2 and maxcut >= m/2 always). So in the uniform case the LP is strictly below the integer odd-cycle answer: 2 vs 3.
- Approximation transfer: if a bipartition S has internal(S) <= alpha*F, then max-internal(S) <= alpha*F <= 2 alpha g. So any alpha-approximation for weighted F (Min UnCut) gives a 2*alpha-approximation for g.

## 1. Q1: judicious partitions (g) in the literature

- Bollobas, Scott, "Judicious partitions of graphs", Combinatorica 19 (1999). [NA original; statement via survey FT]. Theorem 3.3 (sharp, general graphs): every graph with m edges has a bipartition with e(V1,V2) >= m/2 + sqrt(m/8 + 1/64) - 1/8 and e(Vj) <= m/4 + sqrt(m/32 + 1/256) - 1/16 for j=1,2. Extremal: complete graphs of odd order.
- Bollobas, Scott, "Problems and results on judicious partitions", Random Structures & Algorithms 21 (2002). [NA original; survey FT]. Survey: Bollobas, Scott, "Judicious partitions and related problems", https://people.maths.ox.ac.uk/scott/Papers/bccsurvey.pdf [FT]. Contents used: g(m) <= m/3 (equality only K3); g(m) <= m/4 + c m^{4/5} (BS 1993, Period. Math. Hungar.); Porter 1992: g(m) <= m/4 + sqrt(2m); Porter-Yang 2001: m/4 + sqrt(8m/9); Conjecture 3.7 on g(m) < m/4 + sqrt(m/32) - c m^{1/4} for infinitely many m; Problem 3.10 (open): is there a polynomial algorithm finding a partition with max{e(V1),e(V2)} <= g(m)? Survey states no NP-hardness result for judicious partitioning.
- Alon, Bollobas, Krivelevich, Sudakov, "Maximum cuts and judicious partitions in graphs without short cycles", JCTB 88 (2003) 329-346. [FT] (https://web.math.princeton.edu/~nalon/PDFS/abks5.pdf). Key points:
  - Thm 1.1: girth >= r (r>=4) implies f(G) >= m/2 + c m^{r/(r+1)}. Triangle-free (r=4): f >= m/2 + c m^{4/5}.
  - Thm 1.2: for every m there is a girth-5 graph with f <= m/2 + c' m^{5/6} (Erdos-Renyi polarity graphs, pruned). Explicit family, constant unspecified.
  - Thm 1.3 (the bridge to g): if f(G) = m/2 + delta and delta <= m/30, then g(G) <= m/4 - delta/2 + 10 delta^2/m + 3 sqrt(m). (Constant text garbled in extraction; verify.) Remark: if sqrt(m) << delta = o(m), then g = m/4 - (1-o(1)) delta/2.
  - Thm 1.4: if delta >= m/30 then g <= m/4 - m/100.
  - Random graphs G(n,1/2): f >= m/2 + Omega(m^{3/4}), stated tight.
  - Conj 4.1: girth-r tightness known only for r = 4, 5.
  - No NP-hardness discussion.
- Alon, "Bipartite subgraphs", Combinatorica 16 (1996) 301-311. [FT] (https://web.math.princeton.edu/~nalon/PDFS/bipartite3.pdf). Thm 1.2: triangle-free e > 1 edges implies f(G) >= e/2 + c0 e^{4/5}. Tightness: for every e there is a triangle-free G with f(G) <= e/2 + C0 e^{4/5}. Shearer had e^{3/4}. Remark: computing f is NP-complete (MAX-CUT). Randomized poly-time algorithms exist, derandomizable.
- Erdos problem #581 page, https://www.erdosproblems.com/581 [ABS]. Restates f(m) = m/2 + Theta(m^{4/5}) (Alon 1996). No 3/5, 2/3 or golden-ratio constant appears there.
- Gutin, Nielsen, Yeo, Zhou, "Judicious partitions in edge-weighted graphs with bounded maximum weighted degree", arXiv 2507.05827 (2025). [FT]. This is the weighted judicious result, the closest match to our g(y). Thm 2.1 (k=2): a bipartition (X,Y) with w(X,Y) >= w/2 + d_w/4 and max{w(G[X]), w(G[Y])} <= w/4 + Delta_w/8, where w = total weight, d_w = 2w/|V| (average weighted degree), Delta_w = max weighted degree. Tight for odd complete graphs. Prop 2.1: Delta_w cannot replace d_w in the cut bound. No triangle-free/odd-cycle results, no complexity.
- Ji, Ma, Yan, Yu, "On problems about judicious bipartitions of graphs", arXiv 1701.07162. [ABS] (only abstract and theorem excerpts). Thm 1.4: for real lambda >= 1 there is a bipartition with e(V1)^lambda + e(V2)^lambda minimized to within C(t,2)^lambda + C(t+1,2)^lambda, t = sqrt(m/2 + 1/16) - 1/4. Answers Scott's question for lambda=2, implies max{e(V1),e(V2)} bound. Bollobas-Scott parity-bisection conjecture: counterexamples (complete 3-partite).
- Erdos, Faudree, Pach, Spencer, "How to make a graph bipartite", JCTB 45 (1988) 86-98. [NA] (only via pages below). Triangle-free graphs on n vertices can be made bipartite by deleting (1/18 - eps) n^2 edges. Conjectured n^2/25 (blow-up of C5 shows sharpness) [ABS: open-problem pages].
- Sudakov, "Making a K4-free graph bipartite" (PDF, no journal listed) [FT]. K4-free on n vertices: deletion of at most n^2/9 edges suffices, extremal complete 3-partite. Context only (vertex-count bounds, not m-based).

Verdict Q1: g is exactly the Bollobas-Scott judicious bipartition functional (unweighted), and the weighted version is studied by Gutin et al. (2025). Triangle-free results for g are asymptotic: g = m/4 - (1-o(1)) delta/2 by ABKS Thm 1.3 when f = m/2 + delta. No explicit small triangle-free examples found.

## 2. A6: B_graph vs cheapest odd cycle

Argument (checked by hand):
- Let H be triangle-free with m edges and maxcut(H) < 3m/5. Then F(H) > 2m/5, so g(H) >= F/2 > m/5.
- Weights b: b = 1 on E(H), b = M >= 5 on all other pairs. Then B_graph <= b.1_H / g(H) = m/g(H) < 5.
- Every odd cycle of K_K costs >= 5: pure-H odd cycles have length >= 5 (triangle-free), others contain an M-edge. So the cheapest odd cycle is >= 5 (= 5 if H has a C5). Hence B_graph < cheapest odd cycle. A6 holds.
- Existence of H: Alon 1996 gives, for every m, a triangle-free H with maxcut <= m/2 + C0 m^{4/5}. So the chain works for m > (10 C0)^5. C0 is not explicit.
- Explicit alternative: ABKS Thm 1.2 girth-5 polarity family, maxcut <= m/2 + c' m^{5/6}. Also explicit, constants unspecified.
- Stronger, any odd girth L: Dembo, Montanari, Sen, "Extremal cuts of sparse random graphs", Ann. Probab. 2017 (arXiv 1503.03923) [ABS]: for random gamma-regular graphs maxcut/n = gamma/4 + P* sqrt(gamma/4) + o(sqrt gamma) w.h.p., P* ~ 0.7632. So maxcut/m -> 1/2 as gamma -> infinity, with leading correction P*/sqrt(gamma). Deleting one edge from each of the O_p(1) short odd cycles gives odd girth >= L with the same asymptotic ratio. Then with b = 1 on E(H) and b = M >= L off E(H), B_graph <= m/g(H) -> about 4 < L. So the B_graph vs cheapest-odd-cycle gap can be made as large as wanted (gamma non-explicit; the o_gamma error is unspecified). This last step is my argument, not stated in the sources.
- Small examples: C5 gives m/g = 5, not < 5. Petersen: F = 3 so g <= 3 and m/g >= 5, so not a counterexample. In general, a counterexample needs maxcut < 4m/5 at least (since g <= F), and g close to F. No small explicit example found; I did not search computationally (no code allowed).

Verdict A6: the inequality is true asymptotically using known constructions, and the LP form is even stronger (B_bip <= m/F ~ 2). The specific statement B_graph < cheapest odd cycle is not found in the literature (it is an elementary consequence). No small explicit example found.

## 3. Q2: the LP min{b.y : y(E) - y(delta(S)) >= 1 for all S}

- Barahona, Grotschel, Mahjoub, "Facets of the bipartite subgraph polytope", Math. Oper. Res. 10(2) (1985) 340-358. [ABS] (https://pubsonline.informs.org/doi/10.1287/moor.10.2.340). Shows odd-order complete subgraphs and odd bicycle wheels induce facets of PB(G) (the convex hull of incidence vectors of bipartite subgraphs). Note: the DOI 10.1287/moor.10.2.359 is Gerards 1985 (odd bicycle wheel inequalities), not BGM.
- Gerards, "Testing the odd bicycle wheel inequalities for the bipartite subgraph polytope", MOR 10(2) (1985) 359-360. [ABS]. Polynomial-time separation for odd bicycle wheel inequalities.
- Barahona, Mahjoub, "On the cut polytope", Math. Programming 36 (1986). [FT] (https://www.lamsade.dauphine.fr/~mahjoub/Papers/BM86.pdf). Defines PC(G) (cut polytope) and PB(G) with PC contained in PB. Gives clique-type, odd-cycle-type and wheel-type facets for PC. Does not discuss frustration or edge bipartization by name.
- Fiorini, Hardy, Reed, Vetta, "Approximate min-max relations for odd cycles in planar graphs", Math. Programming. [FT] (https://samuel.fiorini.web.ulb.be/papers/odd_final.pdf). Covering LP (min sum y_v with each odd cycle hit) vs packing LP. Goemans-Williamson covering LP integrality gap <= 9/4 (conjectured 3/2 for vertex version; text rendering ambiguous). NP-hardness of min odd cycle transversal and max odd cycle packing, even planar. Edge version: Kral-Voss factor 2.
- Agarwal, Charikar, Makarychev, Makarychev, "O(sqrt(log n)) approximation algorithms for Min UnCut, Min 2CNF deletion, and directed cut problems", STOC 2005. [NA] (cited in CKKRS [FT] and in Agarwal-Alon-Charikar [FT]). Min UnCut = edge bipartization = our F(y), weighted version.

Observations (mine):
- Our LP = fractional relaxation whose integer solutions are shortest odd cycles. Uniform case: LP value 2 vs integer 3 (see section 0).
- The odd-cycle covering LP (min b.y s.t. y(C) >= 1 for odd C) is a different LP. K4 with y = 1/3 on all edges satisfies every triangle constraint (value 1), but F = 2/3 < 1, so it is not feasible for our LP. So the feasible sets are not nested in the obvious way. I did not determine the containment in either direction.
- I did not find the LP by name in the literature. Ingredients (cut/bipartite subgraph polytopes, max-cut separation, odd-cycle relaxations) are known.

Verdict Q2: partially known (ingredients: BM86, BGM85, odd-cycle relaxations; ACMM05 for approximation). The specific LP and the uniform-case values 2 vs 3 appear to be new to me and not found in the sources.

## 4. Q3: judicious partitions / bipartization in majority-vote or jury contexts

- Searches (4 targeted queries): no result applying judicious partitions, odd-cycle covers or bipartization to bound majority-vote or jury error. Results found but not read: Kaniovski, "Condorcet jury theorem: dependent case"-type papers and "Monotonicity in Condorcet's Jury Theorem with dependent voters" (Springer, Social Choice and Welfare 2007) [NA]. arXiv 2607.14288 "Majority Correctness in Social Networks" returned an empty PDF [NA].
- The judicious-partition and bipartization literature does not mention voting.

Verdict Q3: not found.

## 5. Q4: complexity and approximation

- Max-cut NP-complete: stated in Alon 1996 [FT]; Karp 1972 [NA, from memory].
- Max-cut hardness: Khot, Kindler, Mossel, O'Donnell (paper PDF at cornell) [ABS]. Under UGC + Majority Is Stablest: NP-hard to approximate beyond Goemans-Williamson alpha_GW ~ 0.878567. Under UGC alone: beyond 3/4 + 1/(2 pi) ~ 0.909 (as stated in the abstract/theorems extracted).
- Judicious bipartition (unweighted, decision with budgets k1, k2): Lokshtanov, Saurabh, Sharma, Zehavi, "Balanced judicious bipartition is fixed-parameter tractable", arXiv 1710.05491 / FSTTCS 2017 (LIPIcs 93) [FT via PDF]. JB is FPT in k = k1 + k2 (time 2^{O(k)} n^{O(1)}, exact O-constants garbled). The paper gives no NP-hardness proof for JB. Its hardness remark is for the balanced variant only.
- NP-hardness of computing g (unweighted or weighted): NOT found in any source checked. Plausible but unverified. Survey Problem 3.10 remains open for poly-time algorithms.
- Min UnCut O(sqrt(log n)) (ACMM05) gives an O(sqrt(log n)) approximation for g via section 0.
- Min UnCut hardness: CKKRS (https://www.wisdom.weizmann.ac.il/~robi/papers/CKKRS-CutHardness-CompComp06.pdf) [FT] states UGC-hardness for Min-2CNF-equiv deletion (any constant factor; NP-hard within c sqrt(log log n) under a stronger UG variant). Min UnCut is a special case of that problem family (only "not-equal" constraints), so the statement does not automatically transfer to Min UnCut. I did not confirm a Min UnCut-specific statement.
- Computing B_bip exactly: separation is max-cut. Weighted B_bip for general b is an LP with NP-hard separation (poly-time only with approximate separation or via the LP relaxation's structure). Not resolved.

Verdict Q4: judicious bipartition: FPT in budget; NP-hardness not found; Problem 3.10 open. Min UnCut: O(sqrt(log n)) approx; UGC hardness for the related Min-2CNF deletion.

## 6. Other sources checked (context only)

- Erdos problem page (triangle-free n^2/25, open): https://fanchung.ucsd.edu/erdosprb/erdos/newproblems/K3FreeToBipartite.html [ABS]; https://openproblemgarden.org/op/odd_cycle_transversal_in_triangle_free_graphs [ABS]. Both are edge-deletion, vertex-count bounds, and do not apply to m-based bounds.
- Poljak-Tuza survey (DIMACS 1995): not fetched [NA].
- Edwards 1973 bound: cited in survey and Gutin et al.; not fetched [NA].
- Shearer's e/2 + Omega(e^{3/4}) for triangle-free: cited in Alon 1996 [FT].

## 7. Verdict per claim

- A3 (weighted majority, Markov with g): elementary, verified by hand. Not found as a named result in any source. Verdict: not found (proof is short and self-contained). Framing as "bipartition-based bound on majority error" appears new within the sources I checked.
- A6 (B_graph < cheapest odd cycle): known ingredients (Alon 1996, ABKS 2003, DMS 2017); the B_graph statement follows elementarily. Asymptotic truth is solid. Small explicit counterexample: not found. Verdict: partially known (ingredients), statement not found in the literature.
- A8 (F/2 <= g <= F, and the sandwich): elementary identities; max-cut / judicious equivalence not found stated this way. Verdict: elementary, not found by name.
- LP for B_bip (Q2): partially known (BM86, BGM85, Gerards85, ACMM05 approximation). Exact LP and the uniform values 2 and 3 not found. Verdict: partially known.

## 8. Must-cite list

1. Bollobas, Scott, Combinatorica 19 (1999), Judicious partitions of graphs (Thm 3.3 sharp g bound).
2. Bollobas, Scott, RSA 21 (2002), Problems and results on judicious partitions (survey; Problem 3.10).
3. Alon, Bollobas, Krivelevich, Sudakov, JCTB 88 (2003) 329-346 (Thm 1.2, 1.3, Conj 4.1; f <-> g bridge).
4. Alon, Combinatorica 16 (1996) 301-311, Bipartite subgraphs (triangle-free f = m/2 + Theta(m^{4/5})).
5. Gutin, Nielsen, Yeo, Zhou, arXiv 2507.05827, Judicious partitions in edge-weighted graphs (weighted g bound, Delta_w).
6. Lokshtanov, Saurabh, Sharma, Zehavi, FSTTCS 2017 (arXiv 1710.05491), Balanced judicious bipartition is FPT.
7. Barahona, Grotschel, Mahjoub, MOR 10 (1985) 340-358, Facets of the bipartite subgraph polytope.
8. Barahona, Mahjoub, Math. Programming 36 (1986), On the cut polytope.
9. Agarwal, Charikar, Makarychev, Makarychev, STOC 2005, O(sqrt(log n)) for Min UnCut.
10. Khot, Kindler, Mossel, O'Donnell (Max-cut hardness, alpha_GW under UGC + MIS).
11. Dembo, Montanari, Sen, Ann. Probab. 2017, Extremal cuts of sparse random graphs (maxcut/m -> 1/2 for large gamma).
12. Erdos, Faudree, Pach, Spencer, JCTB 45 (1988), How to make a graph bipartite (context; vertex-count).

## 9. Open or unresolved

- NP-hardness of judicious bipartition (weighted or unweighted).
- Min UnCut-specific UGC hardness.
- Containment between the bipartition LP (F >= 1) and the odd-cycle covering LP.
- A small explicit triangle-free (or odd-girth-5) graph H with maxcut(H) < 4m/5 and g(H) > m/5 (would give a concrete counterexample to B_graph = cheapest odd cycle).
- Constants in ABKS Thm 1.3 and Alon's C0 (extraction garbled; verify in originals).
