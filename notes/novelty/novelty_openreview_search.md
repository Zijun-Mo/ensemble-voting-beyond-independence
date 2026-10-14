# OpenReview 检索（2026-10-13，经公开 API api2.openreview.net/notes/search）

说明：内置浏览器不可用（对话未连接电脑），改用 OpenReview 公开检索接口。每个检索只读到前约 100k 字符（约前 15–30 条结果）。NeurIPS 2026 主会投稿尚未公开，无法检索。

## 检索词
1. majority vote dependence
2. weighted majority vote correlated classifiers
3. worst-case ensemble accuracy dependence
4. LLM judges correlated errors aggregation
5. distributionally robust aggregation weak supervision
6. Condorcet jury theorem correlated voters
7. ensemble disagreement unlabeled bound majority vote error
8. minimax robust weights voting unknown correlation pairwise
9. partial identification ensemble accuracy unlabeled agreement sharp bounds

## 结论
没有检索到任何论文在“已知个体准确率 + 两两依赖信息（上界）下，求最坏情况最优的投票规则/权重”这一问题上与我们重叠；也没有论文出现奇环、二部图、共错图等结构。

## 值得在论文中引用或区分的结果
- An & Dasgupta, "When and How Often is Weighted Majority Vote Optimal Under Log Loss?", UAI 2026。只依赖准确率的加权在测度零集合外都不最优（对数损失下）；不涉及依赖下的最坏情况。可作动机引用。
- Agrawal, Pukdee, Balcan, Ravikumar, "Learning from weak labelers as constraints", ICLR 2025。用每个标注器的误差上界作约束、不假设条件独立；但不做投票规则/权重的最坏情况优化，不用两两信息。
- CARE (ICML 2026) 与其前身 "From Many Voices to One"（NeurIPS 2025 workshops）：用潜因子/MRF 建模评审相关性，参数化模型下的可识别性保证，非分布无关。
- "Nine Judges, Two Effective Votes"（NeurIPS 2026 JUDGe workshop）：相关错误的实证动机。
- "Robust AI Evaluation through Maximal Lotteries"（ICLR 2026 workshop AFAA）：歧义集上的最坏情况偏好聚合，不涉及投票者依赖。
- Karge & Rudolph, "The More the Worst-Case-Merrier" (KR 2022)；Karge et al., "To Lead or to be Led" (AAMAS 2024)：陪审团定理推广，已在先前笔记中记录。
- Theisen et al. NeurIPS 2023；Kim et al. NeurIPS 2024；Feofanov et al. JMLR 2024：基于分歧的多数投票误差界，非紧确、不优化规则。

## 局限
- 每个检索只看了前几十条；API 按相关度排序，长尾结果未逐条查看。
- NeurIPS 2026 主会投稿未公开；ICLR 2027 投稿只出现少量。
