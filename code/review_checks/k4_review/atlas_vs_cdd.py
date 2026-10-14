import json, sys, itertools
from fractions import Fraction as F
import cdd, cdd.gmp as G
atlas=json.load(open(sys.argv[1])); E=list(itertools.combinations(range(4),2))
X=[[(m>>i)&1 for i in range(4)] for m in range(16)]
A=[[1,*x,*(x[i]*x[j] for i,j in E)] for x in X]+[[0]*5+[int(k==j) for k in range(6)] for j in range(6)]
ok=True; tot=0
for wstr,verts in atlas.items():
    w=list(map(int,wstr.split(','))); f=[int(sum(w[i]*(2*x[i]-1) for i in range(4))>0) for x in X]+[0]*6
    g=G.copy_generators(G.polyhedron_from_matrix(G.matrix_from_array([[-f[r]]+A[r] for r in range(22)],rep_type=cdd.RepType.INEQUALITY)))
    mine={tuple(F(x) for x in r[1:]) for r in g.array if r[0]==1}
    theirs={tuple(F(n,v[0]) for n in v[1:]) for v in verts}
    ok&= mine==theirs; tot+=len(theirs)
print('atlas JSON == independent cddlib enumeration for all', len(atlas), 'rule orbits:', ok, '; vertices', tot)
