"""Does the report's Lemma 2 cover every legal certificate?  For every S5-orbit of vertices of the 60 per-rule polyhedra
P_F (all legal certificates are convex combos + rays of these), test exactly whether some permutation of the vertex lies in
one of the report's 11 polyhedra (7 monotone Q_{F_w}, 4 'drop-0' polyhedra), built exactly as in the report's code."""
import json, itertools, sys
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U
STATES=U.STATES; MOMENT=U.MOMENT; A0=U.A_BASE_INT
systems=[]
for cw in [(1,0,0,0,0),(1,1,1,0,0),(1,1,1,1,1),(2,1,1,1,0),(2,2,1,1,1),(3,1,1,1,1),(3,2,2,1,1)]:
    ft=np.array([int(2*sum(cw[i]*int(x[i]) for i in range(5))>sum(cw)) for x in STATES])
    systems.append(('mon'+str(cw),A0,np.r_[ft,np.zeros(10,int)]))
for name,fn in [("const_1",lambda x:1),("dict_4",lambda x:int(x[1]==1)),("maj3_4",lambda x:int(x[1]+x[2]+x[3]>=2)),("apex4_4",lambda x:int(2*x[1]+x[2]+x[3]+x[4]>=3))]:
    ft=np.array([fn(x) if x[0]==0 else 0 for x in STATES])
    pr=[MOMENT[m|1]+MOMENT[(m|1)^31] for m in range(32) if (m&1)==0 and fn(STATES[m])==1]
    if name=="const_1": systems.append(('drop0:'+name,A0,np.r_[ft,np.zeros(10,int)]))
    else: systems.append(('drop0:'+name,np.vstack([A0,np.array(pr)]),np.r_[ft,np.zeros(10,int),np.ones(len(pr),int)]))
PERMS=list(itertools.permutations(range(5)))
orbs=[tuple(F(x) for x in o) for o in json.load(open('non_u248_orbits.json'))]
uncovered=[]
for v in orbs:
    hit=False
    for p in PERMS:
        w=U.apply_perm(v,p); wf=np.array([float(x) for x in w])
        for name,A,b in systems:
            if np.all(A@wf-b>=-1e-9):
                if all(sum(int(A[r,k])*w[k] for k in range(16))>=b[r] for r in range(len(b))): hit=True; break
        if hit: break
    if not hit: uncovered.append(v)
print('non-U248 vertex orbits of the per-rule polyhedra:',len(orbs),'  not contained in any of the 11 polyhedra (any permutation):',len(uncovered))
if uncovered:
    v=uncovered[0]; print('example (a,beta,gamma):',[str(x) for x in v])
    beta=v[1:6]; gam=v[6:]; E=U.EDGES
    w=[beta[i]+sum(gam[e] for e,(a,b) in enumerate(E) if i in (a,b))/2 for i in range(5)]
    print('  antisymmetric weights w_i(Q) =',[str(x) for x in w])
    print('  Q on the 32 states:',[str(v[0]+sum(beta[i]*int(STATES[m][i]) for i in range(5))+sum(gam[e]*int(STATES[m][a]*STATES[m][b]) for e,(a,b) in enumerate(E))) for m in range(32)])
json.dump([[str(x) for x in v] for v in uncovered],open('lemma2_uncovered.json','w'))
