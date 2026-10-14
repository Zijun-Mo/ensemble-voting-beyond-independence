"""Review checks for task-1 v2: (1) the report's interval MILP (solve_and_repair_upper) vs. brute force over all 1882
signed threshold rules (interval-marginal worst-case LP per rule); (2) the review's closed-form K=5 certificate
C248 = min_j [a_j + sum_i (beta_ij^+ q_i^+ - beta_ij^- q_i^-) + sum gamma_ij b_ij] over the 248 certificates of U248,
which must satisfy C248 >= C_full (each term is a feasible interval dual for its rule).
Usage: python3 -I check_task1v2.py <dir with report script> <dir with S26 module> <rule_vertices.json>"""
import sys, json, time, itertools, numpy as np
from fractions import Fraction as F
from scipy.optimize import linprog
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[2])
import task1_interval_milp_verified as T
import D2_K5_U248_completeness_proof_and_verifier as U
K=5; E=list(itertools.combinations(range(K),2)); Z=np.array([[(m>>i)&1 for i in range(K)] for m in range(32)])
rules=[]; seen=set()
CAN=[(1,0,0,0,0),(3,1,1,1,1),(2,1,1,1,0),(3,2,2,1,1),(1,1,1,0,0),(2,2,1,1,1),(1,1,1,1,1)]
for w0 in CAN:
    for w in set(itertools.permutations(w0)):
        for sg in itertools.product((1,-1),repeat=5):
            ww=np.array(w)*np.array(sg); f=tuple(int(ww@(2*z-1)>0) for z in Z)
            if f not in seen: seen.add(f); rules.append(np.array(f,float))
assert len(rules)==1882
C=[c.coefficients() for c in U.make_u248()]
def worst(f,ql,qu,b):
    Aub=np.vstack([Z.T,-Z.T,np.array([Z[:,i]*Z[:,j] for i,j in E])]); bub=np.r_[qu,-ql,b]
    r=linprog(-f,A_ub=Aub,b_ub=bub,A_eq=np.ones((1,32)),b_eq=[1],bounds=(0,None),method='highs'); return -r.fun
rng=np.random.default_rng(3); T0=time.time(); n=0; dm=0; d248=[]
while time.time()-T0<200 and n<60:
    k=rng.integers(4,10); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k))
    qq=rng.uniform(0.1,0.4,5); pi=0.6*pi+0.4*np.prod(np.where(Z==1,qq,1-qq),axis=1)
    q=Z.T@pi; e=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in E])
    if np.any(q>=0.5): continue
    eps=F(int(rng.integers(5,40)),1000)
    qr=[F(round(x,4)).limit_denominator(10**4) for x in q]; er=[F(round(x,4)).limit_denominator(10**4) for x in e]
    ql=[max(F(0),x-eps) for x in qr]; qu=[min(F(1,2),x+eps) for x in qr]
    b=[min(qu[i],qu[j],max(F(0),er[t]+eps)) for t,(i,j) in enumerate(E)]
    res=T.solve_and_repair_upper(K,ql,qu,b,time_limit=10)
    up=float(res['verified_upper'])
    brute=min(worst(f,np.array(ql,float),np.array(qu,float),np.array(b,float)) for f in rules)
    c248=min(float(c[0]+sum(max(c[1+i],0)*qu[i]+min(c[1+i],0)*ql[i] for i in range(5))+sum(c[6+t]*b[t] for t in range(10))) for c in C)
    dm=max(dm,abs(up-brute)); d248.append(c248-brute); n+=1
    assert c248>=brute-1e-9
print(f'{n} random interval instances: max |verified MILP upper - brute force over 1882 rules| = {dm:.2e}; '
      f'C248 - C_full: min {min(d248):.2e}, max {max(d248):.4f}, mean {np.mean(d248):.4f}  [{time.time()-T0:.0f}s]')
