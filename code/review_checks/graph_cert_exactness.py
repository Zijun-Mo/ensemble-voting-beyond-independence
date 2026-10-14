"""Is the optimal deterministic risk R* = min{ q_min , B_graph } ?
B_graph = best weighted-graph certificate = dual (eq. disj) restricted to alpha=beta=0 (pure quadratic Q = sum y_e x_i x_j);
the dictator certificate Q = x_i gives q_i.  Both computed by the same MILP with variables fixed."""
import itertools, time, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
T0=time.time(); rng=np.random.default_rng(5)
def states(K): return np.array(list(itertools.product([0,1],repeat=K)))
def solve(K,Z,q,b,pairs,pure):
    P=len(pairs); nc=1+K+P; reps=[k for k in range(len(Z)) if Z[k,0]==0]; nb=len(reps); n=nc+nb
    feat=lambda z: np.r_[1,z,[z[i]*z[j] for i,j in pairs]]
    rows=[];lo=[]
    for t,r in enumerate(reps):
        x=Z[r]; y=1-x
        for v,s,l in ((x,-1,0),(y,1,1)):
            a=np.zeros(n); a[:nc]=feat(v); a[nc+t]=s; rows.append(a); lo.append(l)
        for v in (x,y):
            a=np.zeros(n); a[:nc]=feat(v); rows.append(a); lo.append(0)
    c=np.r_[1,q,[b[p] for p in pairs],np.zeros(nb)]
    lb=np.r_[-np.inf*np.ones(1+K),np.zeros(P),np.zeros(nb)]; ub=np.r_[np.inf*np.ones(nc),np.ones(nb)]
    if pure: lb[:1+K]=0; ub[:1+K]=0
    res=milp(c,constraints=LinearConstraint(np.array(rows),np.array(lo),np.inf),integrality=np.r_[np.zeros(nc),np.ones(nb)],bounds=Bounds(lb,ub))
    return res.fun if res.success else np.inf
for K in [3,4,5,6,7]:
    Z=states(K); pairs=list(itertools.combinations(range(K),2)); n=nR=bad=badR=0; gaps=[]
    while time.time()-T0<45*(K-2) and n<200:
        q=rng.uniform(0.02,0.3,K)
        if rng.random()<0.6: q=q/q.sum()*rng.uniform(0.3,1.0)
        b={p:min(q[p[0]],q[p[1]])*(rng.uniform(0,0.5) if rng.random()<0.75 else 1.0) for p in pairs}
        full=solve(K,Z,q,b,pairs,False)
        if not np.isfinite(full): continue
        cand=min(q.min(),solve(K,Z,q,b,pairs,True)); n+=1
        inR = q.sum()-(K-2)*q.min()<=1
        if cand>full+1e-9: bad+=1; gaps.append(cand-full); badR+=inR
        nR+=inR
    print(f"K={K}: n={n} (in R: {nR}); min{{q_min,B_graph}} > R* in {bad} cases ({badR} in R); max gap={max(gaps) if gaps else 0:.4f}; t={time.time()-T0:.0f}s",flush=True)
