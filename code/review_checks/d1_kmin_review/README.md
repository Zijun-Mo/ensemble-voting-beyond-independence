# 对 [S24]（K_min = 9 报告）的审查脚本（每次运行 < 5 分钟）
- `S24_appendix_Bgraph_Kmin_9_exact_verify.py`：报告附录脚本原样提取，复跑 25 s，输出与报告一致。
- `enum_Q.py 6`：cddlib(GMP) 精确枚举 Q↑(K_6)：62 个顶点、31 个分数顶点、2 个 S_6 轨道 → `Qfrac_orbits_K6.json`。
- `search_Q.py`：K=6 两个轨道的 B 值（均为 1）；K=7、8 的随机 LP 抽样（只抽到分母 3，说明抽样不可靠，已被下面的完整枚举取代）。
- `check_c5k2.py`：C5⊕K2（K=7）与 Lift A/B（K=8）满足全部奇环约束且紧约束秩满（21/21、28/28），确为顶点。
- `b_c5k2.py`：上述三点 B = 1（HiGHS MILP，原值与对偶界一致到 1e-13）。
- `adj_orbits.py n [秒]`：按 S_n 轨道的邻接分解精确枚举 Q↑(K_n) 的顶点轨道（切锥极射线用 cddlib GMP，比值检验用精确整数），可断点续跑。
  K=6：5 个轨道（2 分数），与 enum_Q 一致；K=7：3 秒完成，9 个轨道、6 个分数（状态 `adj_state_K7.json`）；
  K=8：未完成（`adj_state_K8_partial.json`，23 个已发现轨道，15 个未展开；卡在 |S|=4 的整数顶点切锥）。
- `b_k7_orbits.py`：K=7 全部 6 个分数轨道 B = 1。
- `cmp_lifts.py`：K=8 已发现的分母 >3 轨道与 Lift A/B 比对：除 Lift A/B 外另有 5 个满足匹配条件的新轨道（分母 4、5、5、6、7）。
- `K8_new_orbits_found.txt`：K=8 已发现的分母 >3 分数顶点的具体坐标（Lift A/B 及 5 个新轨道）。
- `S25_user_verify_K8_new.py`：[S25] 所附脚本（原样，末尾仅加了保存向量的一行），2 秒，三个向量在 Q↑ 内、紧约束秩 28、B=1。
- `k8_found_orbits_B.py`（需 `user_k8/user_vecs.pkl` 与 `run/adj_state_K8.json` 的相对路径；输出见 `k8_found_orbits_B_output.txt`）：
  用户三个向量的精确可行性与精确秩（sympy）；与审查枚举的轨道比对（三者分别对应审查发现的分母 4 轨道与两个无零边的分母 5 轨道）；
  K=8 已发现的 19 个分数轨道全部 B=1（含分母 6、7）。另发现一个分母 3、零边为单条匹配边、全部正权 ≥1/3、但不是规范 b*(8,1)/3 的轨道，推翻 [S24] 类别 I 的"唯一锁定为规范轨道"。

## 对 [S27]（K_min = 9 完整证明报告）的审查（第十八份审查）
- `S27_bundle/Bgraph_Kmin_9_K8_complete_verify.py`（报告第五节代码，去掉了 markdown 围栏）：24 s，4 项检查全部通过。注意它只展开了 F02–F19，没有展开 F01 与 4 个整数轨道。
- `adj_orbits.py 8 250 skipint`：审查自己的邻接分解（cddlib GMP 求切锥极射线 + 精确比值检验），把除 F01、K44 整数顶点外的 21 个轨道全部展开，没有新轨道。状态 `adj_state_K8.json`。
- `f01_normaliz.py` / `tangent_normaliz.py K44`：用 Normaliz（PyNormaliz 2.24，精确整数）求 F01（56 行）与 K44 整数顶点（928 行）切锥的全部极射线：52,367 条（15.5 s）与 12,503 条（20.5 s）。
- `expand_from_rays.py`：对这两组射线做精确比值检验，得到全部相邻顶点（F01：52,339 条有界边，命中 18 个轨道；K44：12,475 条，命中 21 个轨道），没有新轨道。规范形用报告中的细化规范形，随机抽样与 40,320 置换的完整规范形交叉核对一致。
- `compare_orbits.py`：报告列出的 19 个分数轨道与审查枚举得到的 19 个逐一相同。
- 结论：Q↑(K_8) 恰有 23 个顶点轨道（4 个整数、19 个分数），23 个轨道全部展开，邻接封闭。
- 精确 B 证书（第十九份审查）：`make_B_certs.py`（子代理编写；每个结点的 LP 用 Fractions 单纯形精确求解，浮点只用于挑选划分）生成 `exact_B_certs/F01.json`…`F19.json`。证书是二叉决策树，叶子给出有理乘子。`verify_B_certs.py`（审查方独立编写，纯 Fraction）0.4 s 校验全部 19 份：共 1,256 个叶子，F01 最大（364 个叶子、深度 35），最大分母 285。
