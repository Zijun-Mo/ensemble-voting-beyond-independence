# 依赖不确定下的稳健加权投票

DSAA2011 Task B · Project 2 "Ensemble Voting Beyond Independence" · 四人组队

**研究问题**。集成成员（例如多个 LLM 评审）的错误相关，但相关结构未知。我们只知道各成员的错误率 $`q_i`$，以及成对共同出错概率的上界 $`b_{ij}`$（可以从无标签数据的分歧率估计）。在这种情况下：

- 哪条投票规则的最坏情况错误率最小？
- 这个最小值是多少？
- 能否从有限样本出发，给出可认证的保证？

**目标**。课程报告按 NeurIPS 论文的标准来写，争取课程项目的 3 分加分。暂定题目：*Certifiable Ensemble Voting under Correlated Errors: Exact Minimax Rules, Judicious Partitions, and Certificates from Unlabeled Disagreement*。

> **从哪里读起：**[`选题文档.md`](选题文档.md) 是项目的唯一主文档，包括记号、全部结果、剩余任务、分工和时间表。本 README 只做导航。

---

## 当前进展

理论部分的结论都已确立。剩下的工作是真实数据实验、证明蓝图、文献核对、论文写作，以及可选的 Lean 形式化（见主文档第 4 节）。

| 结果 | 内容 | 状态 |
|---|---|---|
| 1 | 课程 Part I/II：只知准确率、只知两两相关时的界 | 已证（课程要求） |
| 2 | 强对偶；最优规则总是带符号加权投票；$`R_0`$ 内权重非负 | 已证 |
| 3 ★ | 主定理：零截距证书精确，在 $`R_0=\{\sum q_i\le1\}`$ 内 $`R^*=\min\{q_{\min},B_{\rm graph}\}`$；$`B_{\rm graph}`$ 是明智二部划分问题；K=5 时有 27 项闭式 | 已证 |
| 4 ★ | 边界：走出 $`R_0`$ 即失效；奇圈界在 $`K\ge9`$ 时不够（$`K_{\min}=9`$）；只有二部图约束时无法优于单个模型 | 已证 / 计算机辅助（精确） |
| 5 | 精确值是有限个仿射函数的最小值；项数有超指数下界 | 已证 |
| 6 ★ | $`B_{\rm graph}`$ 的计算：精确 MILP、2 倍松弛；分数色数阈值定理，精确计算是 NP 难的 | 已证 |
| 7 | K=3 有 4 项、K=4 有 12 项全参数闭式 | 计算机辅助（精确） |
| 8 ★ | K=5 有 248 项全参数闭式，每项都不可删；附截距分解与两侧显式界 | 计算机辅助（精确） |
| 9 | 确定性规则与随机化规则何时相等、何时不等 | 已证 / 计算机辅助（精确） |
| 10 ★ | 有限样本认证：图认证、置信集认证、K=5 闭式认证、计算间隙、双路径定理 | 已证 |
| 11 ★ | 对称参数（成员错误率与成对上界都相同）下任意 K 的闭式：随机化规则已证；确定性规则的上界闭式只有 $`O(K^2)`$ 项，奇偶 K 结构不同 | 随机化与边界部分已证；确定性 K≤11 为计算机辅助（精确），一般 K 为猜想 |

"计算机辅助（精确）"指证明依赖计算，但全部判断都用有理数完成，并附有可以独立校验的证书。

---

## 目录结构

```
.
├── 选题文档.md                  主文档（第 1 节为完整记号表，附录 A 为来源索引 [S1]–[S31]）
├── 选题文档_归档.md             唯一的归档：被取代的说法、旧编号对照、旧任务、旧术语 → 新术语对照
├── review_counterexample_report.md   [S8] 逐份审查记录（第 1–22 份）
├── sources/                     [S1]–[S7]、[S13]–[S31]：证明文档、推导报告及其证书包
│   ├── S1/S2/S3_*.pdf           课程 Part I/II 的证明、引用式证明、初学者讲解
│   ├── S4–S7, S13–S24_*.tex     推导报告
│   ├── S26, S27_*.md            K=5 闭式与 K_min 的报告
│   ├── S28–S31_*.tex|pdf        K_min=9 简化证明、有限样本认证、248 项闭式的结构解释、对称参数下任意 K 的闭式
│   └── S17/19/21/22/25/26/27/29/31_bundle/   各报告附带的验证脚本与证书数据
├── code/
│   ├── part12_checks/           [S9b] 课程 Part I/II 的数值核验
│   ├── early_experiments/       [S10] K=3 对偶顶点枚举等早期实验
│   ├── review_checks/           [S9] 审查用核验脚本，按主题分子目录（每个脚本运行 < 5 分钟）
│   └── unlabeled_thread/        [S11] 无标签数据线程的代码、结果与图
├── course_part12/               [S1]–[S3] 的 LaTeX 源文件
├── notes/
│   ├── unlabeled_thread_completion.md   [S11] 识别定理与表格数据实验
│   ├── novelty/                 [S12] 新颖性与文献核查笔记（11 份）
│   └── literature/              早期文献调研（社会选择、概率界、无标签界等）及配图
└── history/                     选题文档 v6–v37 的完整备份；early/ 为最早的两份选题文档
```

**编号约定**。主文档用 [S1]–[S31] 引用来源，完整对照见主文档附录 A。

**路径说明**。旧路径只出现在两处历史记录里：审查记录 [S8] 和 `history/` 中的旧版文档。它们写的是整理前的路径（如 `research_notes/依赖不确定下加权投票研究选题/...`）。对照方法：去掉这个前缀，剩下的部分就是本仓库中的相对路径；`novelty_*.md` 现在在 `notes/novelty/` 下。

**术语说明**。部分来源文件仍用旧术语，例如"奇环""锚定双轻""暴露见证"。它们和主文档用语的对照表在 `选题文档_归档.md` 第四节。

---

## 复现

### 环境

- Python ≥ 3.11，依赖见 [`requirements.txt`](requirements.txt)，安装方法：
  ```bash
  pip install -r requirements.txt
  ```
- 精确多面体计算需要 `pycddlib`（GMP 有理数版，3.x API）和 `PyNormaliz`。只有 $`K_{\min}`$ 的顶点枚举用到它们。
- MILP / LP 用 SciPy 自带的 HiGHS。浮点求解器只用来寻找证书；证书本身都用 `fractions.Fraction` 精确校验。

### 运行约定

- 统一用 `python3 -I` 运行。
- 很多脚本按**当前目录**读写数据文件，所以要先 `cd` 到脚本所在目录再运行。
- 需要导入报告附带模块的脚本，用第一个参数传入相应 bundle 目录。
- 每个子目录的 `README.md`（如有）写明了用法和预期输出。

### 示例

两个关键的计算机辅助证明：

```bash
# 结果 4(b)：K=8 时 19 个分数顶点轨道上 B_graph ≥ 1 的有理决策树证书（约 0.5 秒）
python3 -I code/review_checks/d1_kmin_review/verify_B_certs.py \
    sources/S27_bundle/Bgraph_Kmin_9_K8_complete_verify.py \
    code/review_checks/d1_kmin_review/exact_B_certs
# 预期输出：19 certificates verified; orbits without certificate: []

# 结果 8：其余 1,648 个对偶顶点轨道在全部可行参数上都不低于 U248（约 40 秒）
cd code/review_checks/k5_u248_review
python3 -I dominance_full.py ../../../sources/S26_bundle 0 2000
# 预期输出：exact Farkas certificates verified: 1648/1648; failures 0
```

有限样本认证的代码在 `sources/S29_bundle/`：

- `task1_C248_closed_form.py`：K=5 闭式认证，不需要求解器；
- `task1_interval_milp_verified.py`：一般 K 的区间 MILP，附有理数核验与计算间隙。

### 真实数据（[S11]）

`code/unlabeled_thread/q8*.py`、`q11*.py` 用到 6 个 PMLB 表格数据集：adult、magic、spambase、phoneme、twonorm、churn。数据**不随仓库分发**。请从 [PMLB](https://github.com/EpistasisLab/pmlb) 下载 `<名称>.tsv.gz`，放到 `code/unlabeled_thread/data_real/` 下，或者用环境变量 `PMLB_DATA_DIR` 指定目录。

### LaTeX

- 用 `xelatex` 编译。
- `sources/` 中多数推导报告基于 `ctexart`，需要 TeX Live 的 ctex 宏包。
- [S23]、[S28]、[S30]、[S31] 和 `course_part12/` 用 `fontspec` 加 Noto Serif/Sans CJK SC 字体，不依赖 ctex。

---

## 分工与后续工作

- **分工**：见主文档第 6 节。A 负责理论，B 负责计算，C 负责统计，D 负责实验；每人负责一组结果的写作、蓝图与问答。
- **时间表与检查点**：见主文档第 8 节。
- **剩余任务**：见主文档第 4 节。
  1. 真实数据实验
  2. 证明蓝图与证书包
  3. 论文写作
  4. 文献核对与新颖性检索
  5. Lean 形式化（可选）
