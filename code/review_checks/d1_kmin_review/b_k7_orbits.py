"""B_graph on every matching-valid fractional vertex orbit of Q(K_7) from the complete adjacency enumeration (run/adj_state_K7.json)."""
import json,itertools,numpy as np
from fractions import Fraction as F
from scipy.optimize import milp, LinearConstraint, Bounds
n=7;E=list(itertools.combinations(range(n),2));m=len(E)
def bgraph(cost):
    reps=[S for S in range(1<<n) if not S&1]; nz=len(reps); N=m+nz; R=[];lo=[]
    for t,S in enumerate(reps):
        a=np.zeros(N); d=np.zeros(N)
        for k,(x,y) in enumerate(E):
            if (S>>x)&1 and (S>>y)&1: a[k]=1
            if not (S>>x)&1 and not (S>>y)&1: d[k]=1
        a[m+t]=-1; R.append(a); lo.append(0); d[m+t]=1; R.append(d); lo.append(1)
    r=milp(np.r_[cost,np.zeros(nz)],constraints=LinearConstraint(np.array(R),lo,np.inf),integrality=np.r_[np.zeros(m),np.ones(nz)],bounds=Bounds(0,1),options={'mip_rel_gap':1e-12})
    return r.fun, r.mip_dual_bound
for k in json.load(open('run/adj_state_K7.json'))['orbits']:
    v=[F(x) for x in k.split('|')]
    if all(x.denominator==1 for x in v): continue
    z=[E[i] for i,x in enumerate(v) if x==0]; vs=[a for e in z for a in e]
    print(sorted(set(map(str,v)),key=F),'zero',z,'matching' if len(vs)==len(set(vs)) else 'non-matching','B (primal, dual bound) =',bgraph(np.array([float(x) for x in v])))
