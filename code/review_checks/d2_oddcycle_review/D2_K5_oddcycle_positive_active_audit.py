"""Exact audit for the K=5 odd-cycle collapse and positive-intercept activity theorem.
Requires Python 3.10+, sympy. No floating-point arithmetic is used.
The Lehman/MNI implication is a cited theorem, NOT proved by this script.
"""
from itertools import combinations, permutations, product
from fractions import Fraction as F
from random import Random
import sympy as sp

V = range(5)
E = tuple(combinations(V, 2))
edge = {e: k for k, e in enumerate(E)}
tri = [tuple(edge[e] for e in combinations(S, 2)) for S in combinations(V, 3)]
pent = []
for rest in permutations(range(1, 5)):
    if rest[0] > rest[-1]: continue
    cyc = (0,) + rest
    pent.append(tuple(sorted(edge[tuple(sorted((cyc[i], cyc[(i+1)%5])))]
                             for i in range(5))))
cycles = tri + pent
assert len(cycles) == len(set(cycles)) == 22
M = sp.Matrix([[int(k in C) for k in range(10)] for C in tri])
assert M.det() == -48 and M.inv()*sp.ones(10, 1) == sp.ones(10, 1)/3

covers = set()
for mask in range(1 << 10):
    if all(any(mask & (1 << e) for e in C) for C in cycles):
        covers.add(mask)
minimal = {m for m in covers if all((m ^ (1 << e)) not in covers
                                     for e in range(10) if m & (1 << e))}
parts = set()
for mask in range(1, 1 << 5):
    if not(mask & 1) or mask == 31: continue
    parts.add(sum((1 << k) for k, (i, j) in enumerate(E)
                  if bool(mask & (1 << i)) == bool(mask & (1 << j))))
assert minimal == parts and len(parts) == 15

# For any nonnegative s_i, gamma_ij, express Q_a(S)=H(S)-a(|S|-1).
# All a>=0 satisfying quadratic cover constraints form:
#    [0, min(A,T)] union [1,A] (the second interval if A>=1).
# A comprises 26 positivity caps, 5 singleton/complement caps if needed,
# and 10 pair/triple caps; T=(H(V)-1)/4.
subs = list(range(32))
size = [m.bit_count() for m in subs]

def max_a(s, g):
    H = [sum(s[i] for i in V if m & (1 << i)) +
         sum(g[k] for k,(i,j) in enumerate(E)
             if m & (1 << i) and m & (1 << j)) for m in subs]
    caps = [H[m] / (size[m] - 1) for m in subs if size[m] >= 2]
    caps.extend((H[31 ^ (1 << i)]-1)/3
                for i in V if s[i] < 1)
    caps.extend(max(H[m]-1, (H[31 ^ m]-1)/2)
                for m in subs if size[m] == 2)
    A = min(caps)
    T = (H[31]-1)/4
    a = A if A >= 1 else min(A, T)
    return a, A, T, H

def legal(a, H):
    return (a >= 0 and
            all(H[m]-a*(size[m]-1) >= 0 for m in subs) and
            all(max(H[m]-a*(size[m]-1),
                    H[31^m]-a*(size[31^m]-1)) >= 1
                for m in range(16)))

rng = Random(20261013)
checks = 0
for _ in range(120):
    s = [F(rng.randrange(0, 15), 8) for i in V]
    g = [F(rng.randrange(0, 15), 8) for e in E]
    a, A, T, H = max_a(s, g)
    if a < 0: continue     # no feasible intercept for this (s,g)
    assert legal(a, H)
    assert legal(F(0), H)
    # There can be no admissible larger intercept: scan all 52 cap values.
    for step in [F(1,100), F(1,4), F(1)]:
        assert not legal(a+step, H)
    checks += 1
assert checks > 30
print('PASS: 22 odd cycles; 15 minimal integral covers; 1/3 fractional point')
print('PASS: exact positive-intercept endpoint formula on', checks,
      'rational test instances')
print('THEOREM: Lehman MNI implication and cap formula require analytic proofs.')

# Demonstrate the NEW coefficientwise 27-atom convex-hull decomposition.
# This LP audit is numerical; the theorem in the report uses separation.
import numpy as np
from scipy.optimize import linprog
B = np.zeros((15, 27))
B[:5, :5] = np.eye(5)
for j, C in enumerate(cycles):
    for k in C:
        B[5 + k, 5 + j] = 1.0
assert B.shape == (15, 27)
for _ in range(150):
    vec = np.array([rng.randrange(1, 24) / 20 for _ in range(15)])
    hv = [sum(vec[i] for i in V if m & (1 << i)) +
          sum(vec[5 + k] for k,(i,j) in enumerate(E)
              if m & (1 << i) and m & (1 << j)) for m in subs]
    norm = min(max(hv[m], hv[31 ^ m]) for m in range(16))
    assert norm > 0
    vec /= norm             # thus vec is in D_0, in fact on its boundary
    witness = linprog(np.zeros(27), A_ub=B, b_ub=vec,
                      A_eq=[np.ones(27)], b_eq=[1],
                      bounds=(0, None), method='highs')
    assert witness.success
    assert np.min(vec - B @ witness.x) > -1e-8
print('NUMERICAL SUPPLEMENT: 150 normalized zero-intercept certificates')
print('admit 27-atom coefficientwise convex decompositions.')
