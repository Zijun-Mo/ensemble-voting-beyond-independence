"""Deterministic vs randomized minimax (K=5): randomized value = LP min L_Q s.t. Q>=0, Q(x)+Q(1-x)>=1, gamma>=0.
Compared on random parameters with full effective graph K5 and on the 8 remaining graph types."""
import numpy as np, itertools, time, networkx as nx
from scipy.optimize import linprog
from exact_core import Z, FEAT, EDGES, rstar
T0=time.time(); rng=np.random.default_rng(4)
REPS=[m for m in range(32) if not m&1]
A=np.vstack([-FEAT, -(FEAT[REPS]+FEAT[[31-m for m in REPS]])]); bnd=np.r_[np.zeros(32),-np.ones(16)]
def rand_val(q,b):
    r=linprog(np.r_[1,q,b],A_ub=A,b_ub=bnd,bounds=[(None,None)]*6+[(0,None)]*10,method='highs'); return r.fun
stat={}
while time.time()-T0<200:
    k=rng.integers(2,12); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k)); q=Z.T@pi
    if q.min()<.02 or q.max()>=.5 or q.sum()-3*q.min()>1: continue
    J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
    b=J+(cap-J)*rng.uniform(0,.9,10)**2
    H=nx.Graph(); H.add_nodes_from(range(5)); H.add_edges_from(e for k2,e in enumerate(EDGES) if b[k2]<cap[k2]-1e-12)
    key=('K5' if H.number_of_edges()==10 else 'other')+(' R0' if q.sum()<=1 else ' R\\R0')
    d=rstar(q,b)[0]; r=rand_val(q,b); s=stat.setdefault(key,[0,0,0.0]); s[0]+=1
    if d>r+1e-7: s[1]+=1; s[2]=max(s[2],d-r)
for k,v in sorted(stat.items()): print(f'{k}: {v[0]} instances; deterministic > randomized in {v[1]} (max gap {v[2]:.4f})')
