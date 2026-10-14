"""Exp 11: computer-assisted proof of the K=3 capped closed form.
Primal: max P(|E|>=2) s.t. P(E_i)=q_i, P(E_iE_j)<=ebar_ij.  Dual (D): min y0+sum y_i q_i+sum u_ij ebar_ij
s.t. y0 + sum_{i in S} y_i + sum_{ij subset S} u_ij >= 1{|S|>=2}  for all S subset [3];  u>=0.
D is pointed, so the dual optimum is attained at a vertex.  Enumerate all vertices exactly (fractions)."""
import itertools
from fractions import Fraction as F
import numpy as np
import sympy as sp
subsets = [S for r in range(4) for S in itertools.combinations(range(3), r)]
pairs = [(0, 1), (0, 2), (1, 2)]
rows, rhs = [], []
for S in subsets:
    row = [1] + [1 if i in S else 0 for i in range(3)] + [1 if (i in S and j in S) else 0 for (i, j) in pairs]
    rows.append(row); rhs.append(1 if len(S) >= 2 else 0)
for k in range(3):  # u_k >= 0
    row = [0] * 7; row[4 + k] = 1; rows.append(row); rhs.append(0)
A = sp.Matrix(rows); b = sp.Matrix(rhs)
verts = set()
for B in itertools.combinations(range(len(rows)), 7):
    M = A.extract(list(B), list(range(7)))
    if M.det() == 0: continue
    x = M.LUsolve(b.extract(list(B), [0]))
    if all((A.row(r) * x)[0] >= b[r] for r in range(len(rows))):
        verts.add(tuple(x))
q1, q2, q3, e12, e13, e23 = sp.symbols('q1 q2 q3 e12 e13 e23')
terms = sorted({sp.simplify(v[0] + v[1] * q1 + v[2] * q2 + v[3] * q3 + v[4] * e12 + v[5] * e13 + v[6] * e23) for v in verts}, key=str)
print(len(verts), "vertices; objective terms:")
for t in terms: print("  ", t)
