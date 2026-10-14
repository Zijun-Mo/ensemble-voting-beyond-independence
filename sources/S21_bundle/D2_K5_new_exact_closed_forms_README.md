# D2: certified K=5 sparse-graph exact closed forms (2026-10-13)

The relevant model has 5 binary error indicators with means q_i and pairwise error probability upper bounds b_ij. The informative graph contains edge ij iff b_ij < min(q_i,q_j). Assume the ambiguity set of joint distributions is nonempty and 0 <= q_i <= 1/2.

## 1. C5 plus one chord
Graph H is a subset of E = {01, 02, 04, 12, 23, 34}. Then

R_det,5 = T_5 = min { q_0,...,q_4, b_01+b_02+b_12, b_01+b_12+b_23+b_34+b_04 }.

The previously hypothesized extra condition involving q_0+q_1+q_2 is NOT needed.

Run: `python D2_K5_chorded_cycle_complete_exact.py` (~15s)
Files: `D2_K5_chorded_cycle_complete_exact.py` and `D2_K5_chord_34_exact_Farkas.json`.

Proof: exact double-description in 13 homogeneous coordinates reconstructs all 44 irredundant facets of the 54-inequality mixed-dual polyhedron (plus t>=0), from 34 candidate vertices and 24 recession rays. In addition, 34 exact-rational Farkas multipliers prove each vertex is dominated throughout the full feasible moment parameter region by one of the 7 listed candidate costs. Recession rays have nonnegative cost by moment feasibility. The two odd-cycle majority rules and dictators achieve the same upper bound for deterministic voting.

## 2. Three-page book graph
Graph H is a subset of E = {01, 02, 03, 04, 12, 13, 14}. Then

R_det,5 = T_5 = min { q_0,...,q_4, b_01+b_02+b_12, b_01+b_03+b_13, b_01+b_04+b_14 }.

Run: `python D2_K5_book_complete_exact.py` (~10s)
Files: `D2_K5_book_complete_exact.py` and `D2_K5_book_16_exact_Farkas.json`.

Proof: exact double-description in 14 homogeneous coordinates yields 49 irredundant facets and certifies 16 mixed-dual vertices plus 27 recession rays as a complete list. 16 exact-rational Farkas multipliers prove dominance by the 8 displayed affine costs.

## 3. Randomization gap and symmetric mixed exact formula
Even inside the high-precision region R, deterministic and randomized K=5 minimax can differ. For q_i=1/5 and b_ij=1/25,

R_det = 3/25 = 0.12,  T_5 = 1/10 = 0.10.

In general for homogeneous q_i=q and b_ij=b,

T_5(q,b) = min{ q, (5/2)b, 1-5q+10b } (for q<=1/2 and feasible parameters),

whereas previously proven deterministic formula is R_det=min{q,3b,1-5q+10b}.

Run: `python D2_K5_randomization_gap_exact.py`. It checks all 120 possible rank-three active bases of the symmetric mixed dual via exact fractions, obtains exactly six vertices, derives the three active terms for q<=1/2, checks a matching primal witness at q=1/5,b=1/25, and establishes an open heterogeneous parameter box with a strict gap.

These results DO NOT imply randomized and deterministic values agree for arbitrary informative graphs. Remaining unresolved five-node biconnected structures require separate analysis.
