"""Numerical sanity check of [S23] Proposition for K=6: inside the open sets U_C (every odd cycle C) and V_i
(every dictator), the exact deterministic minimax (full no-big-M MILP over all quadratic certificates) equals
the claimed affine piece. Also prints O_K and K+O_K for K=3..9."""
import itertools, numpy as np, time
from math import comb, factorial
from scipy.optimize import milp, LinearConstraint, Bounds
K=6; P=list(itertools.combinations(range(K),2)); Z=np.array([[(m>>i)&1 for i in range(K)] for m in range(1<<K)])
F=np.hstack([np.ones((len(Z),1)),Z,np.array([Z[:,i]*Z[:,j] for i,j in P]).T]); nc=F.shape[1]
REPS=[m for m in range(1<<K) if not m&1]; nb=len(REPS); N=nc+nb; R=[];lo=[]
for t,m in enumerate(REPS):
    r=np.zeros(N); r[:nc]=F[m]; r[nc+t]=-1; R.append(r); lo.append(0)
    r=np.zeros(N); r[:nc]=F[(1<<K)-1-m]; r[nc+t]=1; R.append(r); lo.append(1)
R=np.array(R); lb=np.r_[-np.inf*np.ones(1+K),np.zeros(len(P)),np.zeros(nb)]; ub=np.r_[np.inf*np.ones(nc),np.ones(nb)]
def rdet(q,b): return milp(np.r_[1,q,b,np.zeros(nb)],constraints=LinearConstraint(R,lo,np.inf),integrality=np.r_[np.zeros(nc),np.ones(nb)],bounds=Bounds(lb,ub),options={'mip_rel_gap':1e-10}).fun
cycles=[]
for L in range(3,K+1,2):
    for S in itertools.combinations(range(K),L):
        for perm in itertools.permutations(S[1:]):
            if perm[0]>perm[-1]: continue
            cyc=(S[0],)+perm; cycles.append(frozenset(tuple(sorted((cyc[i],cyc[(i+1)%L]))) for i in range(L)))
cycles=sorted(set(cycles),key=len)
OK=lambda k: sum(comb(k,L)*factorial(L-1)//2 for L in range(3,k+1,2))
print('O_K, K+O_K for K=3..9:',[(k,OK(k),k+OK(k)) for k in range(3,10)])
assert len(cycles)==OK(K)
rng=np.random.default_rng(0); T0=time.time(); bad=0; n=0
for C in cycles:
    q=rng.uniform(3/(4*K),1/K,K); b=np.array([rng.uniform(0,1/(4*K*K)) if e in C else rng.uniform(1/(2*K),5/(8*K)) for e in P])
    v=rdet(q,b); n+=1; bad+= abs(v-sum(b[P.index(e)] for e in C))>1e-8
for i in range(K):
    q=rng.uniform(.90/K,.95/K,K); q[i]=rng.uniform(.40/K,.45/K); b=rng.uniform(.30/K,.39/K,len(P))
    v=rdet(q,b); n+=1; bad+= abs(v-q[i])>1e-8
print(f'K=6: {len(cycles)} odd-cycle open sets + {K} dictator open sets checked ({n} points): mismatches {bad}   [{time.time()-T0:.0f}s]')
