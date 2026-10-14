import numpy as np, time, sys
from exact_core import *
from D2_K5_98_term_closed_form_candidate import make_certificates
C=make_certificates()
COEF=np.array([[float(v) for v in c.coefficients()] for c in C]); LAB=[c.label for c in C]
def U(q,b):
    v=COEF@np.r_[1,q,b]; k=int(np.argmin(v)); return v[k],LAB[k]
T0=time.time(); rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
stats={}; worst=[]
def test(tag,q,b):
    u,lab=U(q,b); r,x=rstar(q,b)
    s=stats.setdefault(tag,[0,0,0,0.0]); s[0]+=1
    if u>r+1e-7:
        s[1]+=1; s[3]=max(s[3],u-r)
        worst.append((u-r,tag,np.round(q,4).tolist(),np.round(b,4).tolist(),u,lab,r,np.round(x,4).tolist()))
    if u<r-1e-7: s[2]+=1   # would be a bug (U is a valid upper bound)
inR=lambda q: q.max()<0.5 and q.min()>0 and q.sum()-3*q.min()<=1+1e-12
# A: old d2_hetero distribution (uniform q in R\R0)
rngA=np.random.default_rng(21); nA=0
while nA<300 and time.time()-T0<80:
    q=rngA.uniform(0.2,0.45,K)
    if not (q.sum()>1 and inR(q)): continue
    b=np.array([min(q[i],q[j])*rngA.uniform(0,0.6) for i,j in EDGES])
    if not feasible(q,b): continue
    nA+=1; test('A_uniform_R\\R0',q,b)
# B: random sparse joint laws -> boundary-rich feasible points
while time.time()-T0<230:
    k=rng.integers(2,12); S=rng.choice(32,k,replace=False); w=rng.dirichlet(np.ones(k)*rng.choice([0.3,1,3]))
    pi=np.zeros(32); pi[S]=w; q=Z.T@pi
    if not inR(q): continue
    J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
    mode=rng.integers(3)
    b=J if mode==0 else J+(cap-J)*rng.uniform(0,1,P)**3*(rng.uniform(size=P)<0.5)
    test('B_sparse_'+('R0' if q.sum()<=1 else 'R\\R0'),q,b)
for k,v in stats.items(): print(f"{k}: n={v[0]}  U98>R* in {v[1]}  (max gap {v[3]:.4g})  U98<R* (bug) {v[2]}")
worst.sort(reverse=True)
for w in worst[:4]: print(w)
print('time',round(time.time()-T0,1))
