"""Independent checks: fan 10-term theorem and K4-subdivision 9-term theorem on R, against the exact
deterministic MILP; counts how often the positive-intercept fan term is the active minimum; exact example check."""
import numpy as np, itertools, time
from scipy.optimize import linprog
from exact_core import Z, FEAT, EDGES, rstar
T0=time.time(); rng=np.random.default_rng(31)
REPS=[m for m in range(32) if not m&1]
A=np.vstack([-FEAT,-(FEAT[REPS]+FEAT[[31-m for m in REPS]])]); bnd=np.r_[np.zeros(32),-np.ones(16)]
def rand_val(q,b): return linprog(np.r_[1,q,b],A_ub=A,b_ub=bnd,bounds=[(None,None)]*6+[(0,None)]*10,method='highs').fun
def bd(b,e): return b[EDGES.index(tuple(sorted(e)))]
FAN=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]; SUB=[(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,4)]
def fan10(q,b):
    t=[*q, bd(b,(0,1))+bd(b,(0,2))+bd(b,(1,2)), bd(b,(0,1))+bd(b,(0,4))+bd(b,(1,4)), bd(b,(0,2))+bd(b,(0,3))+bd(b,(2,3)),
       sum(bd(b,e) for e in [(0,3),(0,4),(1,2),(1,4),(2,3)]), 1-q[0]-q[1]-q[2]+sum(bd(b,e) for e in FAN)]
    return min(t), int(np.argmin(t))
def sub9(q,b):
    t=[*q, sum(bd(b,e) for e in [(0,2),(0,4),(2,4)]), sum(bd(b,e) for e in [(1,2),(1,4),(2,4)]),
       sum(bd(b,e) for e in [(0,2),(0,3),(1,3),(1,4),(2,4)]), sum(bd(b,e) for e in [(0,3),(0,4),(1,2),(1,3),(2,4)])]
    return min(t), int(np.argmin(t))
def sample(H, tri=None):
    kind=rng.integers(5)
    if kind==0: p=rng.uniform(.05,.49,5); pi=np.prod(np.where(Z==1,p,1-p),axis=1)
    elif kind==1: p=np.full(5,rng.uniform(.05,.33)); pi=np.prod(np.where(Z==1,p,1-p),axis=1)
    elif kind==2:
        pi=np.zeros(32)
        for _ in range(rng.integers(2,4)):
            p=rng.uniform(.02,.49,5); pi+=rng.uniform()*np.prod(np.where(Z==1,p,1-p),axis=1)
        pi/=pi.sum()
    elif kind==3:
        k=rng.integers(2,12); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k))
    else:   # targeted: heavy-error triangle {0,1,2} with mostly 'exactly one or all-equal' patterns, light outer vertices
        w=rng.dirichlet(np.ones(6)); pi=np.zeros(32)
        pats=[0b00011,0b00101,0b00110,0b00001,0b00010,0b00100]
        for wi,m in zip(w,pats):
            for ext in range(4):
                x3,x4=ext&1,ext>>1; pi[m|(x3<<3)|(x4<<4)]+=wi*[.49,.17,.17,.17][ext]
        pi=.85*pi+.15*rng.dirichlet(np.ones(32))
    q=Z.T@pi
    if q.max()>=.5 or q.min()<=0 or q.sum()-3*q.min()>1: return None
    J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
    slack=rng.choice([0,0,1])*rng.uniform(0,1,10)**2
    return q, np.where([e in H for e in EDGES], J+(cap-J)*slack, cap)
for name,H,f,deadline in (('fan',FAN,fan10,140),('K4-subdivision',SUB,sub9,250)):
    n=bad=0; act={}; nR0=0
    while time.time()-T0<deadline:
        s=sample(H)
        if s is None: continue
        q,b=s; v,k=f(q,b); n+=1; nR0+=q.sum()<=1; act[k]=act.get(k,0)+1
        if abs(rstar(q,b)[0]-v)>1e-7: bad+=1
    groups={'dictator':sum(act.get(i,0) for i in range(5)),'odd cycles':sum(act.get(i,0) for i in range(5,9))}
    if name=='fan': groups['positive-intercept term']=act.get(9,0)
    print(f'{name}: {n} instances in R ({n-nR0} with sum q>1): exact MILP != closed form in {bad}; active term counts {groups}')
q=np.array([599/2000,199/400,199/400,599/2000,599/2000]); bf={(0,1):797/8000,(0,2):797/8000,(0,3):1/1600,(0,4):1/1600,(1,2):797/8000,(1,4):1/1600,(2,3):1/1600}
b=np.array([bf.get(e,min(q[e[0]],q[e[1]])) for e in EDGES])
print(f'example: MILP deterministic = {rstar(q,b)[0]:.7f} (11/1600 = {11/1600:.7f}); randomized LP = {rand_val(q,b):.7f} (3/640 = {3/640:.7f}); fan10 = {fan10(q,b)[0]:.7f}   [{time.time()-T0:.0f}s]')
