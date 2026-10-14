# 第三轮文献检索：图结构上的成对信息概率界与闭式（Haiku 子代理，2026-10-13）
问题：剩余 216 种标号有效图（8 种结构）或更大图族上，"单变量 + 成对上界"稳健投票是否有现成闭式。
## 子代理结论（A 为已打开来源核实，B 为推断）
- A：Boros, Scozzari, Tardella, Veneziani, "Polynomially Computable Bounds for the Probability of the Union of Events", MOR, doi 10.1287/moor.2014.0657。给定单点与至多 m 阶交的精确概率，得到多项式可算界；假设精确 p_ij，不是上界。
- A：Boros & Lee, "Boole's probability bounding problem, LP aggregations, and nonnegative QPBFs", arXiv 2110.10672。LP 对偶 = 非负二次伪布尔函数锥，与 cut polytope 相连；与我们的证书结构最接近，但未给出精确图类。
- A：Padberg 1989, "The Boolean quadric polytope", Math. Prog. 45:139–172。串并联图上 BQP 的完整线性描述。
- A：Michini, "Tight cycle relaxations for the cut polytope"（optimization-online 2016/06/5483）。转述：CUT(G)=MET(G) ⇔ G 无 K5 子式（Barahona–Mahjoub）；BQP(G)=BQP3(G) ⇔ G 无 K4 子式；BQP(G) ↔ CUT(G+apex)。
- A：Deza & Dutour Sikirić, arXiv 1501.05407。CUT_6 有 368 个面（3 个轨道），完整列表已知。
- A：DRO 方向（Li et al. NeurIPS 2022 树结构矩模糊集；arXiv 2310.20657）与本问题不匹配。
- B：8 种剩余结构中，不含 K4 子式的恰为 house（五环加弦）、gem（扇形）、book K_{1,1,3} 三种（我们已用 networkx 复核结构列表，与之吻合）。
- 底线：存在框架文献（BQP / cut polytope、非负二次伪布尔锥），但没有"成对上界 + 图结构 + 多数投票"的现成闭式；多数文献假设精确 p_ij，而我们只有上界。注意：内层最坏分布问题依赖完整联合分布，相关多面体刻画的是可行域，不直接给出 minimax 值。
Sources: https://pubsonline.informs.org/doi/10.1287/moor.2014.0657 ; https://arxiv.org/pdf/2110.10672v2 ; https://optimization-online.org/wp-content/uploads/2016/06/5483.pdf ; https://arxiv.org/pdf/1506.06702v1 ; https://arxiv.org/pdf/1501.05407 ; https://arxiv.org/pdf/1710.07576
