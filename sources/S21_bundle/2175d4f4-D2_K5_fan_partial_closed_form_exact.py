"""Exact K=5 fan-graph odd-cycle cover certificate and rational example.

Proof checks: all 330 rank-7 potential active bases of the odd-cycle
covering polyhedron; each feasible vertex is a 0-1 monochromatic-edge cut.
This verifies the only graph-specific lemma used in the conditional theorem.

Requires sympy; no floating point in the certification.
"""
import itertools
from fractions import Fraction as F
import sympy as S

E = [(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
C = [ (0,1,4), (0,3,5), (1,2,6), (2,3,4,5,6) ]
A = S.Matrix([[int(e in c) for e in range(7)] for c in C] + [[int(e==k) for e in range(7)] for k in range(7)])
h = S.Matrix([1]*4 + [0]*7)

vertices = set()
for I in itertools.combinations(range(11), 7):
    B = A.extract(I, range(7))
    if B.det(method='domain-ge') == 0:
        continue
    v = B.inv()*h.extract(I, [0])
    if all(t >= 0 for t in A*v-h):
        vertices.add(tuple(v))
cuts = {
    tuple(int(((mask>>i)&1)==((mask>>j)&1)) for i,j in E)
    for mask in range(32)
}
assert len(vertices)==9, len(vertices)
assert vertices.issubset(cuts)
assert all(set(v).issubset({0,1}) for v in vertices)
print('CERTIFIED: all 330 active bases; 9 vertices; each a genuine cut agreement vector')

# Exact feasible heterogeneous instance with R0 violated; cap b on E
q = [F(29,100), F(30,100), F(31,100), F(35,100), F(36,100)]
a = F(4,5)
tri = q[:3]
tri_law = {0: 1-sum(tri)}
for i,p in enumerate(tri): tri_law[1<<i]=p
pi = []
for mask in range(32):
    low = mask & 7
    base = a*tri_law.get(low,F(0))
    for i in (3,4): base *= q[i] if (mask>>i)&1 else 1-q[i]
    indep = (1-a)
    for i,p in enumerate(q): indep *= p if (mask>>i)&1 else 1-p
    pi.append(base+indep)
assert sum(pi)==1
assert all(sum(pi[m] for m in range(32) if (m>>i)&1)==q[i] for i in range(5))
b = {e:sum(pi[m] for m in range(32) if (m>>e[0])&1 and (m>>e[1])&1) for e in E}
t = min([*q,*[sum(b[E[e]] for e in cycle) for cycle in C]])
max_clique_surplus = max(sum(q[i] for i in triangle)-1 for triangle in [(0,1,2),(0,1,4),(0,2,3)])
assert t == F(2699,50000), (t, b)
assert t >= max_clique_surplus
assert sum(q)>1
assert sum(q)-3*min(q)<1
print('q =',list(map(str,q)))
print('b along E =',[str(b[e]) for e in E])
print('t =',t,'=',float(t))
print('triangle surplus maximum =',max_clique_surplus)
print('CERTIFIED: feasible heterogeneous example in R\\R0, exact conditional fan closed-form value t')
