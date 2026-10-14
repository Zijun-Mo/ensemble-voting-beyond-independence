"""Exact (GMP rational) vertex enumeration of the per-rule dual polyhedra
P_F = {(a,beta,gamma): gamma>=0, Q(x)>=F(x) for all 32 x}
for one representative of each of the 60 S5-orbits of signed self-dual threshold rules on 5 voters."""
import itertools, json, time, numpy as np
from fractions import Fraction
import cdd, cdd.gmp as G
from exact_core import Z, FEAT, EDGES
CANON=[(1,0,0,0,0),(3,1,1,1,1),(2,1,1,1,0),(3,2,2,1,1),(1,1,1,0,0),(2,2,1,1,1),(1,1,1,1,1)]
PERMS=list(itertools.permutations(range(5)))
def rule(w):
    f=[]
    for z in Z:
        d=sum(w[i]*(1 if z[i] else -1) for i in range(5)); assert d!=0; f.append(int(d>0))
    return tuple(f)
def perm_rule(f,p):
    out=[0]*32
    for m in range(32):
        img=sum(((m>>i)&1)<<p[i] for i in range(5)); out[img]=f[m]
    return tuple(out)
reps={}; seen=set()
for w0 in CANON:
    for w in sorted(set(itertools.permutations(w0))):
        for sg in itertools.product((1,-1),repeat=5):
            ww=tuple(w[i]*sg[i] for i in range(5)); f=rule(ww)
            if f in seen: continue
            orb={perm_rule(f,p) for p in PERMS}; seen|=orb; reps[f]=(ww,len(orb))
assert len(seen)==1882 and len(reps)==60
T=time.time(); out=[]; tot=0
for f,(w,osz) in reps.items():
    rows=[[-f[m]]+[int(v) for v in FEAT[m]] for m in range(32)]
    rows+=[[0]+[0]*6+[1 if k==e else 0 for k in range(10)] for e in range(10)]
    poly=G.polyhedron_from_matrix(G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY))
    g=G.copy_generators(poly); assert not g.lin_set
    V=[[str(Fraction(x)) for x in r[1:]] for r in g.array if r[0]==1]
    out.append(dict(w=list(w),orbit=osz,rule=list(f),vertices=V)); tot+=len(V)
print('rules',len(out),'total vertices (orbit reps of rules)',tot,'time %.1fs'%(time.time()-T))
json.dump(out,open('rule_vertices.json','w'))
