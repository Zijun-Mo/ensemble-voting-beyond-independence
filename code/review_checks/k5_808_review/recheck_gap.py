import numpy as np, itertools, time
from fractions import Fraction as Fr
from math import comb
from scipy.optimize import linprog
from exact_core import Z, FEAT, EDGES, rstar
REPS=[m for m in range(32) if not m&1]
A=np.vstack([-FEAT, -(FEAT[REPS]+FEAT[[31-m for m in REPS]])]); bnd=np.r_[np.zeros(32),-np.ones(16)]
def rand_val(q,b): return linprog(np.r_[1,q,b],A_ub=A,b_ub=bnd,bounds=[(None,None)]*6+[(0,None)]*10,method='highs').fun
# 1) the user's example, exact checks
q=np.full(5,.2); b=np.full(10,1/25)
print('example q=1/5,b=1/25: deterministic', round(rstar(q,b)[0],6), ' randomized', round(rand_val(q,b),6))
ok=all(Fr(comb(k,2),4)+Fr(comb(5-k,2),4)>=1 for k in range(6)); print('  Q=C(N,2)/4 valid randomized certificate:',ok,'; objective',Fr(1,4)*10*Fr(1,25))
pi={m:(Fr(1,100) if bin(m).count('1') in (2,3) else Fr(1,10) if bin(m).count('1')==1 else Fr(3,10) if m==0 else Fr(0)) for m in range(32)}
qq=[sum(p for m,p in pi.items() if m>>i&1) for i in range(5)]; jj=[sum(p for m,p in pi.items() if m>>i&1 and m>>j&1) for i,j in EDGES]
sym=sum(min(pi[m],pi[31-m]) for m in range(32))
print('  witness: total',sum(pi.values()),' q',set(qq),' pair joints',set(jj),' complementary-symmetric mass',sym)
# 2) why did my sampler miss it?  same sampler vs a dense (mixture-of-product) sampler
rng=np.random.default_rng(4)
def sparse():
    k=rng.integers(2,12); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k)); return pi
def dense():
    pi=np.zeros(32)
    for _ in range(rng.integers(1,4)):
        p=rng.uniform(.02,.45,5); w=rng.uniform()
        pi+=w*np.prod(np.where(Z==1,p,1-p),axis=1)
    return pi/pi.sum()
T0=time.time()
for name,gen in (('sparse (old sampler)',sparse),('dense mixtures of products',dense)):
    n=g=0; mx=0
    while time.time()-T0<(110 if name.startswith('sparse') else 230):
        pi=gen(); q=Z.T@pi
        if q.min()<.02 or q.max()>=.5 or q.sum()-3*q.min()>1: continue
        J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
        b=J+(cap-J)*rng.uniform(0,.9,10)**2
        d=rstar(q,b)[0]; r=rand_val(q,b); n+=1
        if d>r+1e-7: g+=1; mx=max(mx,d-r)
    print(f'{name}: {n} instances in R, deterministic > randomized in {g} (max gap {mx:.4f})')
