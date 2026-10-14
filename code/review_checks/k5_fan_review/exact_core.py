import numpy as np, itertools
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
K=5; EDGES=list(itertools.combinations(range(K),2)); P=len(EDGES)
Z=np.array([[(m>>i)&1 for i in range(K)] for m in range(32)])
FEAT=np.array([np.r_[1,z,[z[i]*z[j] for i,j in EDGES]] for z in Z],float)  # 32x16
REPS=[m for m in range(32) if not (m&1)]
def rstar(q,b):
    """exact R*_det via no-big-M MILP; returns value and (a,beta,gamma)"""
    nc=1+K+P; nb=16; n=nc+nb; rows=[];lo=[]
    for t,m in enumerate(REPS):
        a=np.zeros(n); a[:nc]=FEAT[m]; a[nc+t]=-1; rows.append(a); lo.append(0)
        a=np.zeros(n); a[:nc]=FEAT[31-m]; a[nc+t]=1; rows.append(a); lo.append(1)
    c=np.r_[1,q,b,np.zeros(nb)]
    lb=np.r_[-np.inf*np.ones(1+K),np.zeros(P),np.zeros(nb)]; ub=np.r_[np.inf*np.ones(nc),np.ones(nb)]
    r=milp(c,constraints=LinearConstraint(np.array(rows),np.array(lo),np.inf),integrality=np.r_[np.zeros(nc),np.ones(nb)],bounds=Bounds(lb,ub),options={'mip_rel_gap':1e-9})
    return r.fun, r.x[:nc]
def feasible(q,b):
    A_eq=np.vstack([np.ones(32)]+[Z[:,i] for i in range(K)]); A_ub=np.array([Z[:,i]*Z[:,j] for i,j in EDGES])
    r=linprog(np.zeros(32),A_ub=A_ub,b_ub=b,A_eq=A_eq,b_eq=np.r_[1,q],bounds=(0,None),method='highs'); return r.status==0
def fixed_risk(f,q,b):
    A_eq=np.vstack([np.ones(32)]+[Z[:,i] for i in range(K)]); A_ub=np.array([Z[:,i]*Z[:,j] for i,j in EDGES])
    r=linprog(-f,A_ub=A_ub,b_ub=b,A_eq=A_eq,b_eq=np.r_[1,q],bounds=(0,None),method='highs'); return -r.fun
