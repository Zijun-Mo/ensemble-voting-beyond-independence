# D2: An exact ten-term deterministic minimax theorem for the five-node fan

## Model and assumptions

Five error indicators `X_i∈{0,1}` (indices 0,...,4) have known marginals `q_i` and upper bounds `E[X_i X_j]≤b_ij`. Assume the moment ambiguity set is nonempty, `0≤q_i≤1/2`, and the high-accuracy region `R: sum_i q_i−3 min_i q_i≤1`. The only informative pairwise bounds belong to the fan edge set

`H={01,02,03,04,12,14,23}`.

All other pairwise upper bounds equal their vacuous values `min(q_i,q_j)`. The minimax rule is deterministic and self-dual.

## Exact closed form (new theorem)

`R*_det = min{`

- `q_0,...,q_4;`
- `b_01+b_02+b_12;`
- `b_01+b_04+b_14;`
- `b_02+b_03+b_23;`
- `b_03+b_04+b_12+b_14+b_23;`
- `1−q_0−q_1−q_2+b_01+b_02+b_03+b_04+b_12+b_14+b_23`.

`}`

The first nine entries are dictators and odd-cycle majority rules. The tenth is realized by a weighted majority with integer weights `(2,1,1,2,2)` and any self-dual tie-break.

The tenth entry arises from the nonnegative quadratic majorant

`Q_star(x)=1−x_0−x_1−x_2 + Σ_{e∈H} x_i x_j`.

For the first three vertices, `1−x_0−x_1−x_2+x_0x_1+x_0x_2+x_1x_2` is the indicator that their bits are all equal. The remaining four edges `03,04,14,23` form a path `1−4−0−3−2`; if these edges were all bichromatic, the triangle bits would all be equal. Thus every complementary pair `(x,1−x)` has `max(Q_star(x),Q_star(1−x))≥1`.

## Rigorous completeness proof (computer-assisted finite proof)

1. The mixed/self-dual quadratic dual uses 13 real coefficients, with 32 nonnegativity constraints, 16 complement-pair sum constraints, and 7 nonnegative quadratic coefficients.
2. Homogenize by adding a scale coordinate. An exact integer double-description algorithm verifies that 43 candidate points (including one redundant midpoint) plus 28 recession rays generate the *entire* homogenized cone. The generated cone has 49 facets, each among the 56 defining inequalities. Hadamard's determinant bound certifies that ranks determined modulo two primes are exact over Q.
3. Every deterministic admissible majorant belongs to one of two halfspace branches, `Q(00111)≥1` or `Q(11000)≥1` (bitmask indices 7 and 24).
4. Every vertex of the mixed polyhedron intersected with one of these halfspaces is either an old feasible vertex, a crossing of two old vertices, or a crossing of an old vertex with a recession ray. Exhaustively generating *all* such crossings, including nonedges as harmless extras, produces 239 distinct admissible candidates per branch.
5. For every one of these 478 rational candidates, a sparse exact Farkas witness establishes that its expected quadratic upper bound is at least the ten-term formula *simultaneously over the entire feasible parameter region R*. The check uses only integer/rational arithmetic.
6. Nonnegative recession directions cannot improve the objective when the moment ambiguity set is nonempty. Together with the ten explicit self-dual upper bounds, this proves equality.

## Exact strict-R example

`q=(599/2000,199/400,199/400,599/2000,599/2000)`.

Along edge order `01,02,03,04,12,14,23`, caps are

`b=(797/8000,797/8000,1/1600,1/1600,797/8000,1/1600,1/1600)`.

- `sum q−3 min q = 199/200 <1`, with all `q_i<1/2`.
- Old nine-term expression: `807/8000 = 0.100875`.
- Exact deterministic minimax: `11/1600 = 0.006875`.
- Exact mixed minimax: `3/640 = 0.0046875`.
- Strict deterministic–mixed gap: `7/3200`.

The example script gives a rational joint law, five complementary-pair submasses, and an exact mixed dual upper bound. Strict deterministic optimality follows from the general verified ten-term theorem.

## Verification

Run from the directory containing the files:

```bash
python D2_K5_fan_10term_exact_verify.py
python D2_K5_fan_10term_rational_example.py
```

The theorem verifier requires Python standard library only; no SciPy, solver or external web services. Files:

- `D2_K5_fan_10term_exact_verify.py`: complete integer double-description and exact Farkas checks.
- `D2_K5_fan_10term_generators.json`: 43 homogenized points and 28 recession rays.
- `D2_K5_fan_10term_Farkas_witnesses.json`: all 478 sparse rational dominance witnesses.
- `D2_K5_fan_10term_rational_example.py`: exact counterexample to nine terms and exact mixed gap.

**Scope:** This proof solves the fan effective graph on `R` only, not arbitrary dense K=5 graphs or the general parameter domain outside `R`.
