import time, itertools
from fractions import Fraction as F
import cdd, cdd.gmp as G
from exact_core import FEAT
T=time.time(); REPS=[m for m in range(32) if not m&1]
rows=[[0]+[int(v) for v in FEAT[m]] for m in range(32)]+[[-1]+[int(a+b) for a,b in zip(FEAT[m],FEAT[31-m])] for m in REPS]+[[0]*7+[int(k==e) for k in range(10)] for e in range(10)]
g=G.copy_generators(G.polyhedron_from_matrix(G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY)))
V=[tuple(F(x) for x in r[1:]) for r in g.array if r[0]==1]
E=list(itertools.combinations(range(5),2))
def canon(c):
    best=None
    for p in itertools.permutations(range(5)):
        beta=[0]*5; gam=[0]*10
        for i in range(5): beta[p[i]]=c[1+i]
        for e,(i,j) in enumerate(E): gam[E.index(tuple(sorted((p[i],p[j]))))]=c[6+e]
        t=tuple([c[0]]+beta+gam); best=t if best is None or t<best else best
    return best
orb={canon(v) for v in V}
print(f'randomized K=5 dual polyhedron: {len(V)} vertices, {len(orb)} S5-orbits, intercept>0: {sum(v[0]>0 for v in V)}  [{time.time()-T:.0f}s]')
