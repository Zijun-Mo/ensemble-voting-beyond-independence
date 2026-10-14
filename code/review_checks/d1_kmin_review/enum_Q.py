"""Exact vertex enumeration of the unbounded odd-cycle cover polyhedron Q(K_n) = {u>=0 : u(C)>=1 for all odd cycles C}
(cddlib, GMP rationals); fractional vertices grouped into S_n orbits. Usage: python3 enum_Q.py n"""
import sys, itertools, time, json
from fractions import Fraction as F
import cdd, cdd.gmp as G
n=int(sys.argv[1]); E=list(itertools.combinations(range(n),2)); idx={e:k for k,e in enumerate(E)}
cyc=set()
for L in range(3,n+1,2):
    for S in itertools.combinations(range(n),L):
        for p in itertools.permutations(S[1:]):
            if p[0]>p[-1]: continue
            c=(S[0],)+p; cyc.add(frozenset(idx[tuple(sorted((c[i],c[(i+1)%L])))] for i in range(L)))
rows=[[-1]+[1 if k in C else 0 for k in range(len(E))] for C in cyc]+[[0]+[1 if k==j else 0 for k in range(len(E))] for j in range(len(E))]
T=time.time()
g=G.copy_generators(G.polyhedron_from_matrix(G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY)))
V=[tuple(F(x) for x in r[1:]) for r in g.array if r[0]==1]
frac=[v for v in V if any(x.denominator>1 for x in v)]
perms=list(itertools.permutations(range(n)))
pmaps=[[idx[tuple(sorted((p[i],p[j])))] for i,j in E] for p in perms]
def canon(v):
    best=None
    for pm in pmaps:
        w=[None]*len(E)
        for k,t in enumerate(pm): w[t]=v[k]
        w=tuple(w)
        if best is None or w<best: best=w
    return best
orbs=sorted({canon(v) for v in frac})
print(f'K={n}: {len(cyc)} odd cycles; {len(V)} vertices, {len(frac)} fractional, {len(orbs)} fractional S_{n}-orbits   [{time.time()-T:.0f}s]')
json.dump([[str(x) for x in o] for o in orbs],open(f'Qfrac_orbits_K{n}.json','w'))
