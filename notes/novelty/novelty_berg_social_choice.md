# Novelty check: worst-case voting under pairwise co-error bounds

Research question (as given): K binary voters, accuracies p_i > 1/2, known upper bounds P(Z_i=0, Z_j=0) <= ebar_ij on pairwise co-error, joint distribution otherwise adversarial. Find the rule maximizing worst-case accuracy. Claims: (1) best deterministic rule >= max(p_max, 1 - c*), c* = min over odd cycles of summed ebar on cycle edges, achieved by unweighted majority on that cycle; (2) if the low-co-error graph is bipartite, no rule beats the best single voter; (3) conjecture: the bound is exactly optimal.

Access key: [FT] = full text (or full manuscript) read via fetch, with the content summarized by the fetch tool, so quotes are not verbatim-checked. [ABS] = abstract or record only. [NA] = not accessible.

Searches run: IDEAS/RePEc, Springer, arXiv, Cambridge, UCI site, MPRA, SEP (jury-theorems), and web search for citing work. Google Scholar "cited by" was not available. I found no paper that cites Berg (1994) in a way that states its content.

---

## 1. Papers named in the brief

**Berg (1994), "Evaluation of some weighted majority decision rules under dependent voting", Math. Soc. Sci. 28(2):71-83.** [NA]
- IDEAS record exists, with no abstract. Full text is ScienceDirect-only.
- No keywords. Not cited in the visible SEP bibliography.
- Unknown: dependence model, whether it is worst-case, odd-cycle or bipartite results. The title ("evaluation of some ... rules") suggests a comparison of specific rules under a specific dependence model rather than a maxmin optimization. That is inference from the title only.
- This is the largest unverified item. If Berg 1994 does maxmin over pairwise-constrained dependence for weighted majority, the verdict moves toward "full overlap".

**Shapley & Grofman (1984), "Optimizing group judgmental accuracy in the presence of interdependencies", Public Choice 43(3):329-343.** [FT] (UCI-hosted PDF, read via fetch)
- Independence case: log-odds weights w_i = log(p_i/(1-p_i)) are optimal (their Theorem II).
- Dependence case: no general optimal-weights result. Their Section 3 uses specific joint distributions:
  - Example 1: three voters, exactly one wrong in each state, others right with prob r. They note that for this structure, a rule that opposes the majority on a unanimous vote beats majority, and that the weighted-majority family is not optimal.
  - Example 2: nine voters in three independent districts of three.
  - Theorem III (symmetric distributions): weighted voting is optimal iff p(k) > p(N-k) for all k > N/2. No explicit weights.
  - Section 4: partition voters into blocs with dependence within blocs and independence across them.
- No maxmin/worst-case analysis, no cycle or bipartite structure, no pairwise co-error bounds.
- Overlap with (1)/(2)/(3): low. Relevant to our framing: under dependence the optimal rule can be non-weighted, which supports including non-weighted rules in the rule class.

**Berg (1993), "Condorcet's jury theorem revisited", European J. Political Economy 9(3):437-446.** [ABS, title and metadata only]
- Note a correction: the brief cites "Berg 1993" as a dependency paper. The IDEAS record for this 1993 EJPE paper has no abstract. Its content is unverified.
- The SEP bibliography lists a different item: Berg (1993), "Condorcet's Jury Theorem, Dependency among Jurors", Social Choice and Welfare 10(1). I could not find that record on IDEAS or in web search, so it remains unverified.
- SEP attributes to Berg (1993) the point that "more independence improves majority decisions" and that crowds can still be "wise" under dependence. Overlap: low to medium, since this is a general-dependence comparison, not a worst-case analysis.

**Ladha (1995), "Information pooling through majority-rule voting: Condorcet's jury theorem with correlated votes", J. Econ. Behav. & Org. 26(3):353-372.** [ABS, metadata only]
- No abstract on IDEAS. SEP attributes to Ladha (1992, 1995) a model based on pairwise correlations and a finding that more independence helps. Overlap: low to medium. It uses a specific correlation model, not worst-case.
- Ladha (1992), AJPS, "The Condorcet Jury Theorem, Free Speech, and Correlated Votes": cited in SEP. I did not verify the AJPS record directly. [NA]

**Boland, Proschan & Tong (1989), "Modelling dependence in simple and indirect majority systems", J. Appl. Prob. 26(1):81-88.** [ABS, from a summary of the record]
- Two dependence models are introduced (not named in the record). Reliability-theory framing. Overlap: low, since it fixes specific models.
- SEP separately cites Boland (1989), "Majority Systems and the Condorcet Jury Theorem", JRSS-D 38(3). Not checked.

**Peleg & Zamir (2012), "Extending the Condorcet Jury Theorem to a general dependent jury", Social Choice and Welfare 39(1):91-125.** [NA, no abstract retrieved]
- The title suggests a general dependence condition for CJT asymptotics. Not a worst-case analysis as far as the title indicates. Needs the full text to rule out overlap.

---

## 2. Other social-choice and jury-theorem papers found

**Zaigraev & Kaniovski (2009), "Bounds on the competence of a homogeneous jury" (manuscript).** [FT, manuscript PDF, not the published version]
- Setting: homogeneous jury, exchangeable Bernoulli votes, marginal p, unknown pairwise and higher-order correlations, joint distribution via Bahadur's representation.
- Main results: worst-case min and max of Condorcet probability over all admissible correlation structures, computed as linear programs. Min competence attained at x2 = (2p-1)/(2p). Max = min{2np/(n+1), 1}. Simple majority only; no weighted rules; no pairwise co-error bounds; no odd cycles.
- Overlap: partial. This is the closest methodological match: a worst-case LP over all joint distributions with dependence constraints for majority. It differs in using correlation coefficients (homogeneous, all orders), not co-error upper bounds, and in not optimizing the rule.

**Kaniovski & Zaigraev (2009/2011), "Optimal jury design for homogeneous juries with correlated votes", Theory and Decision 71(4):439-459.** [FT, manuscript PDF]
- Correlation-robust quota alpha = p(n-1)+1 exists, under Bahadur-truncated pairwise correlation. Shows simple majority can be worse than a single juror at low p and high correlation.
- Not worst-case over all dependence structures. Homogeneous. General k-of-n rules, no weights. No pairwise co-error bounds. Overlap: low-partial.

**Kaniovski (2009), "An invariance result for homogeneous juries with correlated votes", Math. Soc. Sci. 57(2):213-222.** [ABS, paraphrase]
- Distributions that keep Condorcet's probability invariant to second-order correlation. Homogeneous juries. Overlap: low.

**Kuncheva, Whitaker, Shipp & Duin (2003), "Limits on the majority vote accuracy in classifier fusion", Pattern Analysis & Applications 6:22-31.** [ABS, from the record and a paraphrase]
- Derives upper and lower bounds on majority-vote accuracy as a function of individual accuracy p, number of classifiers L, and pairwise Yule's Q. The bounds appear to range over all joint distributions consistent with these quantities (extremal "success"/"failure" patterns).
- Caveat: a related ICPR 2000 preprint (lkICPR00.pdf, [ABS]) described the analysis as exhaustive enumeration for three classifiers and small D, not closed-form bounds. The two descriptions differ, and I could not reconcile them without the journal full text. The 2003 journal version is more relevant.
- Uses Q (a correlation statistic), not co-error upper bounds. Majority only. No weighted rules, odd cycles, or bipartite structure. Overlap: partial, in the ML community. Mainly relevant as prior art for "worst case over joint distributions given pairwise dependence, majority vote".

**Pivato (2016/17), "Epistemic democracy with correlated voters", MPRA 69546 (published as J. Math. Econ. per SEP, not re-checked).** [ABS, summary of record]
- Condition: average correlation shrinks with population size. Rule is "a suitably chosen voting rule", unspecified. No worst-case, no pairwise co-error bounds, no odd cycles. Overlap: low.

**Dietrich & Spiekermann (2013), "Independent opinions? On the causal foundations of belief formation and jury theorems", Mind 122(487):655-685.** [ABS]
- Four notions of opinion independence derived from a causal network. No worst-case analysis; no voting-rule optimization. Overlap: none to low.

**Berend & Paroush (1998), "When is Condorcet's Jury Theorem valid?", Social Choice and Welfare 15(4):481-488.** [ABS, from record]
- Necessary and sufficient conditions for CJT under independence with heterogeneous competence. Independence case only. Overlap: none.

**Berend & Sapir (2007), "Monotonicity in Condorcet's Jury Theorem with dependent voters", Social Choice and Welfare 28(3):507-528.** [NA, no abstract]
- Title suggests monotonicity of correctness probability in jury size under dependence. Not about worst-case rules. Overlap: low, unverified.

**Sapir (1998), "The Optimality of the Expert and Majority Rules Under Exponentially Distributed Competence", Theory and Decision 45(1):19-36.** [NA, no abstract]
- Title suggests optimal rules under random competence, not dependence. Overlap: low, unverified.

**Berg (1997), "Indirect voting systems: Banzhaf numbers, majority functions and collective competence", EJPE 13(3):557-573.** [NA, no abstract]
- Title suggests indirect (two-tier) systems. Overlap: low, unverified.

**Braha & de Aguiar (arXiv:2607.14288, July 2026), "Majority Correctness in Social Networks".** [ABS, fetched summary]
- Zealot-contrarian voter dynamics on networks. No worst-case over dependence, no weighted rules, no odd cycles. Overlap: none.

**Gilbert, Schmid, Schnell, Svoboda & Yeo (arXiv:2608.16738, Aug 2026), "Liquid democracy under vote correlation".** [ABS, fetched summary]
- Delegation under signal-dependent competence, with correlation through a shared public signal. Not worst-case over pairwise dependence. Overlap: none to low.

**Aeeneh, Zlatanov & Yu (arXiv:2309.09564, 2023), "New Bounds on the Accuracy of Majority Voting for Multi-Class Classification".** [ABS]
- Independent voters only (i.i.d. and non-i.i.d.). No pairwise dependence, no worst case over joint distributions, no weighted voting. Overlap: none.

---

## 3. Probability-bound literature (mathematical neighbors, not social choice)

These bound P(union) or P(at least k of n events) given partial information. Our c* construction is a union-bound argument on the wrong events, so these are the most important neighbors for claim (1) and for the tightness conjecture (3).

- **Chen (2003), "A lower bound using Hamilton-type circuits and its applications", J. Appl. Prob. 40(4):1121-1132, DOI 10.1239/jap/1067436104.** [ABS] Degree-two lower bound on P(union) using pairwise joint probabilities as edge weights on a Hamilton-type circuit. Structurally close to an odd-cycle bound. Tightness not established in the abstract. Overlap with (1): medium-structural; not about majority voting.
- **Survey: arXiv:1710.07576 (union-probability bounds using partial information).** [FT, fetched summary] Covers Boole, Bonferroni, Kounias, Hunter (spanning-tree upper bound), Dawson-Sankoff, de Caen, KAT, Gallot-Kounias. Dawson-Sankoff is optimal among bounds that depend only on sum P(A_i) and sum P(A_i A_j). Gallot-Kounias uses the full pairwise matrix. The survey contains no cycle or circuit bound.
- **Ramachandra & Natarajan (arXiv:2006.00516, 2020/2022), "Tight Probability Bounds with Pairwise Independence".** [ABS] Tight upper bound on P(at least one of n events) under pairwise independence, with ratio to Boole at most 4/3, and tight at-least-k bounds in special cases. Assumes pairwise independence, not upper bounds on co-error, so it does not overlap directly. But it shows that worst-case, tight probability bounds from pairwise information are an active topic.
- **Not checked, but likely essential:** the correlation polytope and cut polytope literature (Boole 1862; Pitowsky 1989/1991, "Correlation polytopes: their geometry and complexity", Math. Programming; Deza & Laurent 1997, "Geometry of Cuts and Metrics"). From memory, the feasible set of pairwise joint probabilities for binary variables is the correlation polytope, whose facets include odd-cycle inequalities. If so, the odd-cycle structure in (1)/(3) is already the known obstruction to pairwise consistency, and a tightness proof would need to be framed against that literature. The searches I ran returned only quantum-foundations papers in this area. These references are from memory and I did not verify them in this session.

---

## 4. Checks on the brief's citations

- Berg (1993): the brief's "Berg 1993" may refer to the EJPE paper "Condorcet's jury theorem revisited" (IDEAS record found), or to the SCW 10(1) paper cited by SEP. Both remain unverified in content.
- Boland (1989): two different Boland 1989 items appear (JRSS-D; J. Appl. Prob. with Proschan & Tong). Not checked individually.
- Shapley-Grofman: the Springer page shows no optimal-weights result under dependence. The PDF (read via fetch) shows Theorem III (symmetric distributions) and examples. No general weight formula under dependence.

---

## 5. Verdict

**Verdict: PARTIAL overlap.**

Justification: No social-choice paper I could access poses the exact question: a maxmin rule over unknown joint distributions constrained only by pairwise co-error upper bounds, with any rule allowed. The closest social-choice work fixes a dependence model (Ladha; Boland-Proschan-Tong; Berg 1993; Pivato; Dietrich-Spiekermann), or does worst-case over all joint distributions but only for homogeneous juries with correlation coefficients (Zaigraev-Kaniovski LP bounds; Kuncheva et al. 2003 with pairwise Q). Shapley-Grofman (FT) shows that under dependence the optimal rule can be non-weighted, which supports our rule class, but does no worst-case analysis. The probability-bound literature (Chen 2003; Hunter/Kounias/Dawson-Sankoff; Ramachandra-Natarajan) supplies worst-case bounds from pairwise information, so the core union-bound argument behind claim (1) is standard in spirit: majority wrong implies some adjacent pair on an odd cycle is wrong, because an independent set of a cycle of odd length L has size at most (L-1)/2. The novelty is therefore likely concentrated in claim (2) and especially the tightness conjecture (3), which I could not check against the correlation-polytope/Boole literature. The main open risk is Berg (1994), which I could not access: its title suggests an evaluation under a specific dependence model rather than a maxmin, but this is unverified. If its full text does worst-case over pairwise-constrained dependence with weighted majority, the verdict should move toward full overlap.

---

## 6. Recommended next steps

1. Obtain Berg (1994) full text (library access to ScienceDirect, or email the author). This is the single most important check.
2. Read Shapley-Grofman Section 3 directly from the PDF to confirm the Example 1 and Theorem III claims (fetched summary only).
3. Check the correlation-polytope / cut-polytope literature (Pitowsky 1991; Deza & Laurent 1997) for odd-cycle facets and tightness over pairwise-constrained joint distributions. This is the main threat to conjecture (3).
4. Check Chen (2003) full text for whether the Hamilton-circuit bound is tight, and whether it is stated for at-least-k events.
5. Retrieve the Peleg & Zamir (2012) and Berend & Sapir (2007) abstracts through the DOI pages (dx.doi.org) to close the remaining unverified social-choice items.
6. Check the Kuncheva et al. (2003) journal full text to settle the discrepancy with the ICPR 2000 description.

## 更新（2026-10-13，用户提供 ScienceDirect 摘要页截图）
Berg, S. (1994). Evaluation of some weighted majority decision rules under dependent voting. Mathematical Social Sciences 28(2):71–83. DOI 10.1016/0165-4896(94)00752-7. [ABS]
- 摘要：研究相关投票下一类特殊的加权多数规则（如两种专家规则与简单多数）；引入“条件能力概率依赖于一位有影响力成员的意见与能力”的模型（意见领袖模型）；结论为某一相关水平使集体能力最低，负相关可能提高效率。关键词：Correlated votes; Weighted majority rules; Simple majority rule; Expert rules。
- 判断：特定参数化依赖模型 + 固定的少数规则的评估，无最坏情况/max–min、无两两上界约束、无奇环/二部图结构。与我们的问题重叠：低。作为相关工作引用并说明区别即可。
- 参考文献中有 Gradstein et al. (1986) "Performance evaluation of some special classes of weighted majority rules", Math. Soc. Sci.（同一谱系）。被引 23 次，未逐一查看。
