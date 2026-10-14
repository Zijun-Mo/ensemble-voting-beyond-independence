"""Check D2's heterogeneous K=5 conjecture  R*_det = min{q_min, c*, B5}  inside R.
Upper bound on R* is certified independently: take the MILP's optimal certificate, build the
threshold rule, and evaluate its exact worst case with a separate LP."""
import numpy as np, itertools, time
from scipy.optimize import linprog
exec(open('graph_cert_exactness.py').read().split("for K in [3,4,5,6,7]:")[0])
T0=time.time(); K=5; Z=states(K); pairs=list(itertools.combinations(range(K),2))
tri=[c for c in itertools.combinations(range(K),3)]
def cstar(b):
    best=min(b[(a,c)]+b[(a,d)]+b[(c,d)] for a,c,d in tri)
    for rest in itertools.permutations(range(1,5)):
        C=(0,)+rest
        if C[1]>C[-1]: continue
        best=min(best,sum(b[tuple(sorted((C[j],C[(j+1)%5])))] for j in range(5)))
    return best
def risk(f,q,b):
    A_eq=np.vstack([np.ones(len(Z))]+[Z[:,i] for i in range(K)]); A_ub=np.vstack([Z[:,i]*Z[:,j] for i,j in pairs])
    r=linprog(-f,A_ub=A_ub,b_ub=[b[p] for p in pairs],A_eq=A_eq,b_eq=np.r_[1,q],bounds=(0,None),method='highs'); return -r.fun
def best_rule_risk(q,b):
    # re-solve MILP, recover beta,gamma -> threshold rule -> exact LP risk
    from scipy.optimize import milp, LinearConstraint, Bounds
    P=len(pairs); nc=1+K+P; reps=[k for k in range(len(Z)) if Z[k,0]==0]; nb=len(reps); n=nc+nb
    feat=lambda z: np.r_[1,z,[z[i]*z[j] for i,j in pairs]]; rows=[];lo=[]
    for t,r in enumerate(reps):
        x=Z[r]; y=1-x
        a=np.zeros(n); a[:nc]=feat(x); a[nc+t]=-1; rows.append(a); lo.append(0)
        a=np.zeros(n); a[:nc]=feat(y); a[nc+t]=1; rows.append(a); lo.append(1)
    c=np.r_[1,q,[b[p] for p in pairs],np.zeros(nb)]
    lb=np.r_[-np.inf*np.ones(1+K),np.zeros(P),np.zeros(nb)]; ub=np.r_[np.inf*np.ones(nc),np.ones(nb)]
    res=milp(c,constraints=LinearConstraint(np.array(rows),np.array(lo),np.inf),integrality=np.r_[np.zeros(nc),np.ones(nb)],bounds=Bounds(lb,ub))
    beta=res.x[1:1+K]; gam=res.x[1+K:nc]; w=beta.copy()
    for (i,j),g in zip(pairs,gam): w[i]+=g/2; w[j]+=g/2
    best=np.inf
    for s in range(5):
        ww=w+np.random.default_rng(s).normal(0,1e-7,K); f=((2*Z-1)@ww>0).astype(float); best=min(best,risk(f,q,b))
    return res.fun,best
rng=np.random.default_rng(21); n=fails=0; ex=[]
while time.time()-T0<200 and n<300:
    q=rng.uniform(0.2,0.45,K)
    if not (q.sum()>1 and q.sum()-(K-2)*q.min()<=1): continue
    b={p:min(q[p[0]],q[p[1]])*rng.uniform(0,0.6) for p in pairs}
    if not np.isfinite(solve(K,Z,q,b,pairs,False)): continue
    n+=1; cand=min(q.min(),cstar(b),1-q.sum()+sum(b.values()))
    milpv,ruler=best_rule_risk(q,b)
    if ruler<cand-1e-7:
        fails+=1
        if len(ex)<2: ex.append(dict(q=np.round(q,4).tolist(),b={f"{i}{j}":round(v,4) for (i,j),v in b.items()},candidate=round(cand,4),milp=round(milpv,4),certified_rule_risk=round(ruler,4)))
print(f"K=5 instances in R\\R0: {n}; candidate min(q_min,c*,B5) beaten by an LP-certified rule in {fails}")
for e in ex: print(e)
