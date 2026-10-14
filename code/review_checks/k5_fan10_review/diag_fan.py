import numpy as np, itertools, time
from scipy.optimize import linprog
exec(open('check_fan10_sub9.py').read().split("for name,H,f,deadline")[0])
from exact_core import fixed_risk
RULES=None
def truth(q,b):
    global RULES
    if RULES is None:
        R=set()
        for w0 in [(1,0,0,0,0),(3,1,1,1,1),(2,1,1,1,0),(3,2,2,1,1),(1,1,1,0,0),(2,2,1,1,1),(1,1,1,1,1)]:
            for w in set(itertools.permutations(w0)):
                for sg in itertools.product((1,-1),repeat=5): R.add(tuple(int(sum(w[i]*sg[i]*(2*z[i]-1) for i in range(5))>0) for z in Z))
        RULES=[np.array(f,float) for f in R]
    return min(fixed_risk(f,q,b) for f in RULES)
# 1) reproduce the mismatch with the same RNG stream
rng=np.random.default_rng(31); T0=time.time(); found=None; n=0
while time.time()-T0<130 and found is None:
    s=sample(FAN)
    if s is None: continue
    q,b=s; v,k=fan10(q,b); n+=1; r=rstar(q,b)[0]
    if abs(r-v)>1e-7: found=(q,b,v,k,r)
if found:
    q,b,v,k,r=found; t=truth(q,b)
    print(f'mismatch at instance {n}: formula {v:.8f} (term {k}), MILP {r:.8f}, exact (min over 1882 rule LPs) {t:.8f}; q={np.round(q,4)}, sum q={q.sum():.3f}')
# 2) targeted family near the README example: heavy triangle 012, light outer vertices -> positive-intercept term active
def feasible(q,b):
    P=np.array([Z[:,i]*Z[:,j] for i,j in EDGES])
    return linprog(np.zeros(32),A_ub=P,b_ub=b,A_eq=np.vstack([np.ones(32),Z.T]),b_eq=np.r_[1,q],bounds=(0,None),method='highs').status==0
m=act=bad=0; T1=time.time()
while time.time()-T1<140:
    q=np.array([rng.uniform(.27,.33),rng.uniform(.44,.499),rng.uniform(.44,.499),rng.uniform(.27,.33),rng.uniform(.27,.33)])
    if q.sum()-3*q.min()>1: continue
    bf={(0,1):rng.uniform(.05,.15),(0,2):rng.uniform(.05,.15),(1,2):rng.uniform(.05,.15)}
    for e in [(0,3),(0,4),(1,4),(2,3)]: bf[e]=rng.uniform(0,.01)
    b=np.array([bf.get(e,min(q[e[0]],q[e[1]])) for e in EDGES])
    if not feasible(q,b): continue
    v,k=fan10(q,b); m+=1; act+=k==9
    if abs(rstar(q,b)[0]-v)>1e-7 and abs(truth(q,b)-v)>1e-7: bad+=1
print(f'targeted family: {m} feasible instances in R, positive-intercept term active in {act}, mismatches {bad}   [{time.time()-T0:.0f}s]')
