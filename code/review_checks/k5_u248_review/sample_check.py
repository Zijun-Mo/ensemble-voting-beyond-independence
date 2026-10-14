"""Sanity check: exact deterministic minimax (no-big-M MILP, exact_core.rstar) vs U248 on random feasible points of R,
from structured joint distributions (independent, mixtures with sparse/complementary supports); b = p(pi) + random slack."""
import sys, time, itertools, numpy as np
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from exact_core import rstar, Z, EDGES
C=np.array([[float(x) for x in c.coefficients()] for c in U.make_u248()])
rng=np.random.default_rng(7); T=time.time(); n=0; worst=0; nout=0; arg=None
while time.time()-T<240:
    kind=rng.integers(4)
    if kind==0:
        qq=rng.uniform(0.05,0.45,5); pi=np.prod(np.where(Z==1,qq,1-qq),axis=1)
    else:
        k=rng.integers(3,9); S=rng.choice(32,k,replace=False); w=rng.dirichlet(np.ones(k))
        pi=np.zeros(32); pi[S]=w
        if kind==2: qq=rng.uniform(0.1,0.4,5); pi=0.7*pi+0.3*np.prod(np.where(Z==1,qq,1-qq),axis=1)
    q=Z.T@pi
    if np.any(q>=0.5) or np.any(q<=0.01) or q.sum()-3*q.min()>1: continue
    p=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES])
    mq=np.array([min(q[i],q[j]) for i,j in EDGES])
    s=rng.random(10)*(rng.random(10)<0.5)
    b=p+s*(mq-p)
    v,_=rstar(q,b); u=(C[:,0]+C[:,1:6]@q+C[:,6:]@b).min()
    n+=1; nout+=q.sum()>1
    if abs(u-v)>worst: worst=abs(u-v); arg=(q.round(4),b.round(4),u,v)
print(f'{n} feasible points in R ({nout} with sum q>1): max |U248 - R*_det| = {worst:.2e}')
if worst>1e-6: print(arg)
