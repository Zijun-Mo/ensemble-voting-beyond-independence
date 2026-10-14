# Novelty check: Narasimhamurthy (2005/2003) and Mazzetto et al. (AISTATS 2021)

Scope: does either paper already pose or solve "worst-case accuracy of a voting rule given accuracies p_i and pairwise co-error upper bounds ebar_ij" (claims 1 to 3 in the brief)?

## 0. Access legend and what was not reached

- [FT-x] = full text read through WebFetch. WebFetch returns a model-processed extract, not raw text, and quotes are capped at 125 characters. Treat every [FT-x] claim as "extract says", not verbatim.
- [ABS] = abstract or metadata record only.
- [BIB] = bibliographic listing only.

Not reached:
- Narasimhamurthy TPAMI 2005 full text: not accessed. IEEE Xplore pages returned no text (JS-rendered). ACM DL, ResearchGate publication page, PubMed, Semantic Scholar and the Springer pages failed or were rate-limited.
- Narasimhamurthy 2003 "A framework for the analysis of majority voting": not located. The MCS 2003 proceedings volume is LNCS **2709** (verified on Springer), not 2749 as the brief says. The chapter's TOC page 2 was rate-limited (HTTP 429), so I could not confirm the chapter exists in that volume.
- Mazzetto et al. appendix: only the Prop. 5 proof was extracted from the supplementary PDF, and the proof's symbols were partly garbled.
- Mazzetto et al. main paper: read as an extract (see section 2), not raw text. No arXiv version was found.

## 1. Narasimhamurthy, "Theoretical bounds of majority voting performance for a binary classification problem", IEEE TPAMI 27(12):1988-1995, 2005

Identifiers: DOI 10.1109/TPAMI.2005.249. Author: Anand M. Narasimhamurthy (Penn State). Metadata from OpenAlex [ABS]: 17 references, 60 citations.

Abstract, [ABS] via OpenAlex summary. Verbatim quotes that the tool returned:
- "We formulate the majority voting problem as an optimization problem with linear constraints."
- "For a binary classification problem, given the accuracies of the classifiers in the team"

Paraphrase of the rest of the abstract:
- No independence assumption.
- Upper and lower bounds on majority-vote performance are derived from the individual classifier accuracies.
- The objective is nonlinear only for an even number of classifiers when rejection is allowed. Otherwise it is linear, i.e. an LP.
- Two classifier diversity measures are examined against majority-vote performance.

Figure 1 (ResearchGate figure page) [ABS, partial page]. Caption verbatim: "Majority voting accuracy (Odd number of classifiers): Theoretical lower, upper bounds, and independent classifiers for (a) three classifiers, (b) five classifiers, (c) seven classifiers, and (d) nine classifiers with equal accuracy." The page says the bounds are the extremes of majority accuracy over dependence structures consistent with the given accuracies. The constraint set itself is not shown.

Answers to the brief's questions:
- **Information assumed:** accuracies only (the abstract says "given the accuracies"). Pairwise co-error or disagreement is NOT an input. Diversity measures are used only post hoc.
- **LP-based?** Yes, the abstract says so. The even-n case with rejection is nonlinear.
- **Weighted votes or subset choice?** Not in the abstract. The bounds are over the full odd-size team. Abstract-level evidence is absent; not verified in full text. (OpenAlex concept tags include "Weighted voting", but that is an auto-tag and is not evidence.)
- **Closed forms?** Not in the abstract. My own derivation (section 4) gives a closed form for the accuracy-only case with 3 equal classifiers: min majority accuracy = (3p-1)/2 for p in [1/3, 1].

Relevance: this is the accuracy-only analogue. It supports the baseline claim that with no pairwise information the best worst-case is p_max (for claim 1 with c* = infinity). It does not contain any pairwise input, odd cycles, or rule optimality over all Boolean rules.

## 2. Mazzetto, Sam, Park, Upfal, Bach, "Semi-supervised aggregation of dependent weak supervision sources with performance guarantees", AISTATS 2021, PMLR 130:3196-3204

Metadata [ABS] via PMLR and mlanthology (authors match the brief). PDF links verified:
- Main: http://proceedings.mlr.press/v130/mazzetto21a/mazzetto21a.pdf
- Supplementary: http://proceedings.mlr.press/v130/mazzetto21a/mazzetto21a-supp.pdf
- NSF PAR copy: https://par.nsf.gov//servlets/purl/10273522 (main text with references only; appendix not included)
- Code: https://github.com/BatsResearch/mazzetto-aistats21-code
- Undergrad thesis by Dylan Sam (Brown 2021): https://cs.brown.edu/research/pubs/theses/ugrad/2021/sam.dylan.pdf

Abstract (paraphrase of the extract, [FT-x]):
- Theoretical guarantees for learning from weak labelers without assuming independent errors or any particular distribution.
- Efficiently selects small labeler subsets whose majority vote has a provably low error.
- Experiments on 45 image classification tasks.
- About 3 points better than the adversarial method on accurate labelers, and about 15 points better than semi-supervised Dawid-Skene on inaccurate ones.

Setup [FT-x]:
- Known per labeler: error rate eps_i. Known from unlabeled data: pairwise **disagreement** d_ij = P(l_i(x) != l_j(x)). The extract says |eps_i - eps_j| <= d_ij <= eps_i + eps_j.
- Section 4.2 alternative: the full output distribution over {0,1}^n, estimable from unlabeled data.
- Quote (extract): "We make no assumptions on the choice of the labelers in S and possible dependencies between the labelers."

### 2a. Proposition 5 (3-labeler majority), the closed form and its proof

Statement [FT-x, extract]: For I = {i,j,k}, the worst-case error of the majority vote is
  eps_i + eps_j + eps_k - (d_ij + d_jk + d_ik)/2 - 2q,
  with q = 1/2 * max{ eps_i+eps_j-d_ik-d_jk, eps_j+eps_k-d_ij-d_ik, eps_i+eps_k-d_ij-d_jk, 0 }.
Equal case: 3eps - 3d/2 - 2 max{eps-d, 0}. This gives eps + d/2 for d <= eps, and the majority beats eps iff d >= 4eps/3.

Proof [FT-x, appendix extract]: the LP has 8 variables p_abc and 7 equality constraints. Fixing p_000 = x determines the rest. The objective is maximized by the smallest feasible x, and nonnegativity gives the lower bounds on x. The proof says "the objective function is maximized by picking the least x that yields a feasible solution."

**My own check (inclusion-exclusion, to confirm the formula):**
- Let c_ij = P(i and j both wrong) = (eps_i + eps_j - d_ij)/2, since d_ij = P(exactly one wrong).
- P(majority wrong) = P(>=2 wrong) = sum_pairs c_ij - 2t, where t = P(all three wrong).
- Feasibility: t >= 0, t >= c_ij + c_ik - eps_i for each i, and t <= c_ij.
- Worst case: minimize t, giving 2t = max{0, max_i (eps_j + eps_k - d_ij - d_ik)}. This matches Prop. 5 exactly.
- Conclusion: the formula is an elementary derivation and is consistent with the extract. The extract's "3eps - d/2" typo in the summary was an error on my side. The paper's "eps + d/2" is correct.

### 2b. Other propositions [FT-x]
- Prop. 3, verified in the main-paper extract as verbatim text: "Proposition 3. If the only information given to the classification algorithm is the labels of the labelers in ..." Extract says: no function can beat the most accurate labeler, and the worst-case error of the best function equals min eps_i. This is the accuracy-only analogue of claim (1)'s baseline p_max, and it is over all functions, not only majority.
- Props. 1 and 2: constructions showing majority can fail under dependence.
- Section 4.3: a greedy extension that adds two labelers at a time, starting from each labeler.

### 2c. Brief questions
- **Weighted votes:** the main-text extract contains "taking a majority vote weighted by a function of each labeler's accuracy is the optimal way" (Introduction). I believe this refers to the independent case, in the sense of Nitzan-Paroush and Boland. The method itself uses unweighted majority over odd subsets. This sentence was verified only as a raw-text search hit, so the context is unconfirmed.
- **Optimality among all rules:** Prop. 3 covers optimality over all functions only in the accuracy-only case. There is no optimality claim of the form "odd-subset majority is optimal among all rules when pairwise co-errors are known".
- **Odd cycles, bipartite graphs:** absent (term search on the main text). Not in the appendix extract either.
- **Randomized rules:** absent.
- **Linear program:** yes. Conclusion: "We devise a linear program based approach to analytically compute a worst-case error bound of a set of labelers' majority vote."
- **Future work:** absent from the main text. The conclusion has no future-work paragraph. The appendix was not available.
- **Thesis (Sam 2021) [FT-x, extract]:** works only with eps and pairwise disagreement D_ij = P(l_i != l_j). LP over correctness patterns. No separate both-wrong quantity. Future-work sentence: "A future problem of interest is considering alternative methods to use weak supervision sources." Mentions "weighted vote" only via Nitzan-Paroush in related work.

## 3. Related prior art seen (bibliographic only, [BIB]/[ABS])

These are the nearest neighbours for "bounds on majority accuracy with dependence". I did not read their full text.
- Kuncheva, Whitaker, Shipp, Duin, "Limits on the majority vote accuracy in classifier fusion", Pattern Analysis and Applications 6(1):22-31, 2003 [ABS]. Derives upper and lower limits on majority accuracy given individual accuracy, ensemble size, and pairwise dependence (Yule's Q). The strongest neighbour on the dependence-bounds side. It uses pairwise dependence measures, not co-error bounds.
- Ruta and Gabrys, "A theoretical analysis of the limits of majority voting errors for multiple classifier systems", PAA 5(4):333-350, 2002, DOI 10.1007/s100440200030 [ABS]. Claims that a finite set of classifiers can reach 0% majority error even when each has error below 50%, with the exact quote "it is theoretically possible to obtain a 0% majority voting error using a finite number of classifiers." This is relevant to the 5-cycle example (perfect 5-majority), but their construction is distribution-based.
- Boland 1989 (Condorcet jury), Smyth 1996 (PRL, bounds on mean error of multiple experts), Lam and Suen 1997 (SMC), Oh 2003 (PRL, majority accuracy vs dependency), Murua 2002 (PAMI, upper bounds for linear combinations) [BIB].
- Kuncheva, Whitaker, Narasimhamurthy, "A case-study on naive labelling...", Pattern Recognition 2008 [BIB]. Shows Narasimhamurthy's later co-authorship; not on majority bounds.

## 4. Derived observations (my own, not claims from the papers)

1. **Triangle case of claim (1) is elementary.** For any joint distribution, P(majority of 3 wrong) <= sum of the three pairwise both-wrong probabilities (union bound). So the worst-case triangle error is <= c_12+c_23+c_13 <= sum of the ebar values. That is the K=3 case of claim (1).
2. **Odd cycle C_{2m+1}, the general upper bound is also elementary.** Any set S of at least m+1 wrong classifiers must contain a cycle edge, because the independence number of C_{2m+1} is m. So 1[|S|>=m+1] <= number of cycle edges inside S. Taking expectations gives P(majority wrong) <= sum of the cycle co-errors. The 5-cycle example is then a direct instance: zero co-error on cycle edges forces the wrong set to be independent in C5, and so of size at most 2, so 5-majority is never wrong. Mazzetto's framework would not certify the 5-cycle with 3-subsets, since no triangle has full pairwise information.
   **Consequence for the brief:** the upper-bound direction of claim (1) is not the novelty. Novelty has to come from (a) the tightness (equality) conjecture (3), (b) optimality of the threshold rule over all rules (including weighted and randomized), and (c) the minimum over all odd cycles, which needs an LP-duality argument.
3. **Accuracy-only sanity check.** For 3 equal classifiers with accuracy p, the worst-case majority accuracy is (3p-1)/2 for p in [1/3, 1]. Example: p = 0.9 gives 0.85, which is below 0.9. So the best single classifier beats majority in the worst case, consistent with claim (1)'s p_max term and with Mazzetto Prop. 3.
4. **Input difference is the key distinction.** Mazzetto and Narasimhamurthy assume pairwise disagreement or the full distribution, or accuracy only. Our setting has only one-sided co-error bounds ebar_ij (inequalities). Mazzetto's equality constraints do not cover our setting. Under upper bounds, the adversary may choose c_ij below ebar_ij, which changes the worst case.

## 5. Overlap with claims (1) to (3)

- **Claim (1):** partial overlap.
  - Narasimhamurthy: LP worst-case over dependence structures, but accuracy-only input, no pairwise, no cycles, no rule optimality.
  - Mazzetto: closed form for triangle majority with pairwise disagreement (equivalent to exact co-error given accuracies), LP over odd subsets, and Prop. 3 for the accuracy-only baseline. No cycles, no inequality constraints, no optimality over all rules when pairwise info is present.
- **Claim (2)** (bipartite low co-error graph means no rule beats the best single classifier): not covered. Mazzetto Prop. 3 is the accuracy-only analogue and says nothing about co-error constraints. Whether the Prop. 3 adversary respects our edge bounds is not checked.
- **Claim (3)** (conjecture: equality in (1)): not covered by either paper.
- **Example (5-cycle):** not covered by either paper.

**Verdict:**
- Narasimhamurthy 2005: partial (LP worst-case framing and accuracy-only baseline; different input).
- Mazzetto 2021: partial (closest work; pairwise-disagreement LP, closed form for 3-subsets, baseline Prop. 3). Full overlap is not present.
- Claims (2), (3), and the 5-cycle example: none in either paper.

**Must cite:**
- Mazzetto et al. 2021 (closest; cite the triangle closed form and the subset-LP approach, and differentiate on the input type and cycles).
- Narasimhamurthy 2005 (accuracy-only LP bounds; the nearest "worst-case over dependence" framing).
- Kuncheva et al. 2003 (dependence limits on majority accuracy).
- Ruta and Gabrys 2002 (finite-ensemble 0% majority error).
- Nitzan and Paroush 1982 and Boland 1989 (weighted majority optimal under independence).
- Dawid and Skene 1979 and the Dawid-Skene baselines (the independence assumption that Mazzetto relaxes).

## 6. Open items to close before submission
- Obtain the Narasimhamurthy full text (e.g. library access or an author copy) to check whether pairwise terms or weighted votes appear anywhere, and to get the constraint set.
- Confirm whether the 2003 MCS chapter exists in LNCS 2709 (TOC page 2 was rate-limited).
- Read the Mazzetto raw appendix to check the Prop. 5 proof and whether any cycle or weight discussion appears (the extract could have missed it).
- Check the Nitzan-Paroush context in the Mazzetto introduction ("weighted ... is the optimal way").
