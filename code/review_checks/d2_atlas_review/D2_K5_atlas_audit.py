"""Self-contained audit for the finite-atlas and exact K=5 graph-MILP theorems.
Requires sympy, numpy, scipy. The proofs are mathematical; MILP output is
only an implementation check, not an exact rational proof of optimality.
"""
from itertools import combinations
from math import comb
from fractions import Fraction
import sympy as sp
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

K = 5
E = list(combinations(range(K), 2))
S = [tuple((m >> i) & 1 for i in range(K)) for m in range(32)]
ROWS = [[1, *x, *(x[i] * x[j] for i, j in E)] for x in S]
ROWS += [[int(k == j) for j in range(16)] for k in range(6, 16)]
A = sp.Matrix(ROWS)
assert A.shape == (42, 16) and A.rank() == 16

# A fractional positive-intercept vertex for the rule (3,2,2,1,1).
w = (3, 2, 2, 1, 1)
f = sp.Matrix([int(sum(w[i] * (2*x[i] - 1) for i in range(5)) > 0)
               for x in S] + [0] * 10)
z = sp.Matrix([sp.Rational(1,2), sp.Rational(1,2), 0, 0,
               sp.Rational(-1,2), sp.Rational(-1,2),
               0, 0, 0, 0, 0,
               sp.Rational(1,2), sp.Rational(1,2),
               sp.Rational(1,2), sp.Rational(1,2), 1])
slack = A * z - f
assert all(t >= 0 for t in slack)
tight = [i for i,t in enumerate(slack) if t == 0]
assert A[tight, :].rank() == 16
print('EXACT: positive-intercept rational vertex (a=1/2) verified')
print('Uniform vertex numerator/denominator bound:', 2**32)
print('Raw atlas upper bound:', 1882 * comb(42, 16))

# Exact K=5 graph-certificate MILP.
# y_e in [0,1] continuous; z_ij in {0,1} selects the small-side edge.
G, lb = [], []
for i in range(5):
    row = np.zeros(20)
    for k,(u,v) in enumerate(E):
        if i not in (u,v): row[k] = 1
    G.append(row); lb.append(1.0)
for k,(i,j) in enumerate(E):
    row = np.zeros(20); row[k] = 1; row[10+k] = -1
    G.append(row); lb.append(0.0)  # y_ij >= z_ij
    row = np.zeros(20); row[10+k] = 1
    for m,(u,v) in enumerate(E):
        if i not in (u,v) and j not in (u,v): row[m] = 1
    G.append(row); lb.append(1.0)  # y(E(V\\{i,j})) + z_ij >= 1
G = np.asarray(G)
assert G.shape == (25,20)

def graph_value(b):
    result = milp(np.r_[np.asarray(b, dtype=float), np.zeros(10)],
                  integrality=np.r_[np.zeros(10),np.ones(10)],
                  bounds=Bounds(np.zeros(20), np.ones(20)),
                  constraints=LinearConstraint(G, np.array(lb), np.inf),
                  options={'mip_rel_gap':1e-11})
    assert result.success
    y = result.x[:10]
    g = min(max(sum(y[k] for k,(i,j) in enumerate(E) if x[i] and x[j]),
                    sum(y[k] for k,(i,j) in enumerate(E) if not x[i] and not x[j]))
            for x in S)
    assert g > 1-1e-7
    return result.fun

assert abs(graph_value(np.ones(10)) - 3) < 1e-8
print('NUMERICAL: graph MILP, 10 continuous + 10 binary, unit costs -> 3')
print('NOTE: no claim here that arbitrary B_graph equals a shortest-odd-cycle value.')
