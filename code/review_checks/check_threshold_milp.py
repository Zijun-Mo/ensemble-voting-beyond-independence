"""Fast check of the report's Thm 6.1 (optimal self-dual rule = signed linear threshold) and
Thm 7.1 (sum q <= 1 => nonnegative weights suffice).
Exact best value over ALL self-dual rules via ONE MILP on the dual certificate (eq. disj):
  min a + sum beta_i q_i + sum gamma_ij b_ij
  s.t. per complementary pair (x, 1-x): Q(x) >= z, Q(1-x) >= 1-z, Q >= 0, z binary, gamma >= 0.
Then build the threshold rule from the optimal Q (w_i = beta_i + 1/2 sum_j gamma_ij) and
evaluate its worst case by LP; for sum q <= 1 also build the nonnegative-weight rule from Q'."""
import itertools, time, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
T0=time.time(); rng=np.random.default_rng(7)
def states(K): return np.array(list(itertools.product([0,1],repeat=K)))
def risk(f,Z,q,b,pairs):
    A_eq=np.vstack([np.ones(len(Z))]+[Z[:,i] for i in range(Z.shape[1])]); b_eq=np.r_[1,q]
    A_ub=np.vstack([Z[:,i]*Z[:,j] for i,j in pairs]); b_ub=np.array([b[p] for p in pairs])
    r=linprog(-f,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=(0,None),method='highs')
    return None if not r.success else -r.fun
def best_all_rules(K,Z,q,b,pairs):
    P=len(pairs); nc=1+K+P; reps=[k for k in range(len(Z)) if Z[k,0]==0]
    idx={tuple(z):k for k,z in enumerate(Z)}; nb=len(reps); n=nc+nb
    feat=lambda z: np.r_[1,z,[z[i]*z[j] for i,j in pairs]]
    rows=[];lo=[]
    for t,r in enumerate(reps):
        x=Z[r]; y=1-x
        a=np.zeros(n); a[:nc]=feat(x); a[nc+t]=-1; rows.append(a); lo.append(0)      # Q(x)-z>=0
        a=np.zeros(n); a[:nc]=feat(y); a[nc+t]=1;  rows.append(a); lo.append(1)      # Q(y)+z>=1
        a=np.zeros(n); a[:nc]=feat(x); rows.append(a); lo.append(0)
        a=np.zeros(n); a[:nc]=feat(y); rows.append(a); lo.append(0)
    c=np.r_[1,q,[b[p] for p in pairs],np.zeros(nb)]
    lb=np.r_[-np.inf*np.ones(1+K),np.zeros(P),np.zeros(nb)]; ub=np.r_[np.inf*np.ones(nc),np.ones(nb)]
    integ=np.r_[np.zeros(nc),np.ones(nb)]
    res=milp(c,constraints=LinearConstraint(np.array(rows),np.array(lo),np.inf),integrality=integ,bounds=Bounds(lb,ub))
    v=res.x; return res.fun, v[0], v[1:1+K], v[1+K:nc]
def rule_from(beta,gamma,K,pairs,Z,pos=False):
    w=beta.copy()
    for (i,j),g in zip(pairs,gamma): w[i]+=g/2; w[j]+=g/2
    w=w+ (rng.uniform(1e-7,2e-7,K) if pos else rng.normal(0,1e-7,K))   # generic tie-break
    return ((2*Z-1)@w>0).astype(float), w
stats={}
for K in [3,4,5,6,7]:
    Z=states(K); pairs=list(itertools.combinations(range(K),2)); s=dict(n=0,thr_bad=0,nn_n=0,nn_bad=0,neg_needed=0)
    while time.time()-T0 < 40*(K-2) and s['n']<150:
        q=rng.uniform(0.02,0.35,K)
        if rng.random()<0.5: q=q/q.sum()*rng.uniform(0.4,1.0)
        b={p:min(q[p[0]],q[p[1]])*(rng.uniform(0,0.6) if rng.random()<0.7 else 1.0) for p in pairs}
        if risk(np.zeros(len(Z)),Z,q,b,pairs) is None: continue
        val,a,beta,gamma=best_all_rules(K,Z,q,b,pairs); s['n']+=1
        f,w=rule_from(beta,gamma,K,pairs,Z); rv=risk(f,Z,q,b,pairs)
        s['thr_bad']+= abs(rv-val)>1e-7
        if (w<-1e-6).any(): s['neg_needed']+=1
        if q.sum()<=1:
            s['nn_n']+=1; beta2=beta+a   # Q' : alpha'=0, beta'_i=beta_i+a
            f2,w2=rule_from(beta2,gamma,K,pairs,Z,pos=True); assert (w2>=-1e-9).all()
            s['nn_bad']+= abs(risk(f2,Z,q,b,pairs)-val)>1e-7
    stats[K]=s; print(K,s,f"t={time.time()-T0:.0f}s",flush=True)
