"""Independent checks (structured + random instances).
(a) Fan moment-domain theorem (Thm 3.3): closed-form criterion vs LP feasibility of the 5-variable sub-distribution,
    on random (h,r,d) including infeasible ones, plus boundary-perturbed points.
(b) C5+chord and book: exact deterministic MILP and randomized LP vs the closed forms, on structured instances
    (independent products with b = true joint, symmetric q, mixtures, boundary caps) and random ones.
(c) Fan conditional theorem: MILP vs t0 when the triangle-surplus condition holds.
(d) Randomization-gap open box: deterministic MILP vs randomized LP at random points of the box."""
import numpy as np, itertools, time
from scipy.optimize import linprog
from exact_core import Z, FEAT, EDGES, rstar
T0=time.time(); rng=np.random.default_rng(23)
FAN=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
def fan_criterion(h,r,d):
    d=dict(zip(FAN,d)); ok=all(0<=x<=h+1e-12 for x in r) and all(v>=-1e-12 for v in d.values())
    ok&=all(r[i]+r[j]-h<=d[(i,j)]+1e-12 for i,j in [(0,3),(0,4),(1,2),(1,4),(2,3)])
    L01=max(0,r[0]+r[1]-h,r[0]+r[1]+r[4]-h-d[(0,4)]-d[(1,4)]); U01=min(r[0],r[1],d[(0,1)])
    L02=max(0,r[0]+r[2]-h,r[0]+r[2]+r[3]-h-d[(0,3)]-d[(2,3)]); U02=min(r[0],r[2],d[(0,2)])
    ok&= L01<=U01+1e-12 and L02<=U02+1e-12 and L01+L02<=r[0]+d[(1,2)]+1e-12 and U01+U02>=r[0]+r[1]+r[2]-h-d[(1,2)]-1e-12
    return ok
def lp_feasible(h,r,d):
    P=np.array([Z[:,i]*Z[:,j] for i,j in FAN])
    res=linprog(np.zeros(32),A_ub=P,b_ub=d,A_eq=np.vstack([np.ones(32),Z.T]),b_eq=np.r_[h,r],bounds=(0,None),method='highs')
    return res.status==0
agree=dis=nf=0
while time.time()-T0<50:
    h=rng.uniform(.2,1); r=rng.uniform(0,h,5)*rng.choice([1,.5]); d=rng.uniform(0,h/2,7)*rng.choice([1,.3,.1])
    if rng.uniform()<.4:   # move to a near-boundary point: start from a true sub-distribution, shrink caps to the realized joints
        pi=rng.dirichlet(np.ones(32)*rng.choice([.2,1,5]))*h; r=Z.T@pi; J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in FAN]); d=J*rng.choice([1,1-1e-3,1+1e-3],size=7)
    a=fan_criterion(h,r,d); b=lp_feasible(h,r,d); nf+=b
    if a==b: agree+=1
    else: dis+=1
print(f'(a) fan moment-domain criterion vs LP feasibility: {agree+dis} cases ({nf} feasible), disagreements {dis}')
REPS=[m for m in range(32) if not m&1]
A=np.vstack([-FEAT,-(FEAT[REPS]+FEAT[[31-m for m in REPS]])]); bnd=np.r_[np.zeros(32),-np.ones(16)]
def rand_val(q,b): return linprog(np.r_[1,q,b],A_ub=A,b_ub=bnd,bounds=[(None,None)]*6+[(0,None)]*10,method='highs').fun
def structured(H):
    kind=rng.integers(4)
    if kind==0:   # independent product
        p=rng.uniform(.05,.49,5); pi=np.prod(np.where(Z==1,p,1-p),axis=1)
    elif kind==1: # symmetric independent
        p=np.full(5,rng.uniform(.05,.49)); pi=np.prod(np.where(Z==1,p,1-p),axis=1)
    elif kind==2: # mixture of products
        pi=np.zeros(32)
        for _ in range(rng.integers(2,4)):
            p=rng.uniform(.02,.49,5); pi+=rng.uniform()*np.prod(np.where(Z==1,p,1-p),axis=1)
        pi/=pi.sum()
    else:         # sparse
        k=rng.integers(2,12); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k))
    q=Z.T@pi
    if q.max()>=.5: return None
    J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
    slack=rng.choice([0,0,1])*rng.uniform(0,1,10)**2
    b=np.where([e in H for e in EDGES], J+(cap-J)*slack, cap)
    return q,b
def cyc_cost(b,cyc): return sum(b[EDGES.index(tuple(sorted(e)))] for e in cyc)
CH=[(0,1),(0,2),(0,4),(1,2),(2,3),(3,4)]; BOOK=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4)]
forms={'C5+chord':(CH,lambda q,b:min(q.min(),cyc_cost(b,[(0,1),(0,2),(1,2)]),cyc_cost(b,[(0,1),(1,2),(2,3),(3,4),(0,4)]))),
       'book':(BOOK,lambda q,b:min(q.min(),*[cyc_cost(b,[(0,1),(0,k),(1,k)]) for k in (2,3,4)]))}
for name,(H,f) in forms.items():
    n=bd=br=ns=0
    while time.time()-T0<(110 if name=='C5+chord' else 170):
        s=structured(H)
        if s is None: continue
        q,b=s; v=f(q,b); n+=1; ns+=q.sum()>1
        if abs(rstar(q,b)[0]-v)>1e-7: bd+=1
        if abs(rand_val(q,b)-v)>1e-7: br+=1
    print(f'(b) {name}: {n} structured instances ({ns} with sum q>1): deterministic != formula {bd}, randomized != formula {br}')
# (c)
FANE=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
n=bd=0
while time.time()-T0<230:
    s=structured(FANE)
    if s is None: continue
    q,b=s
    t0=min(q.min(),cyc_cost(b,[(0,1),(0,2),(1,2)]),cyc_cost(b,[(0,1),(0,4),(1,4)]),cyc_cost(b,[(0,2),(0,3),(2,3)]),cyc_cost(b,[(0,3),(0,4),(1,2),(1,4),(2,3)]))
    if t0<max(0,q[0]+q[1]+q[2]-1,q[0]+q[1]+q[4]-1,q[0]+q[2]+q[3]-1)-1e-12: continue
    n+=1; bd+= abs(rstar(q,b)[0]-t0)>1e-7
print(f'(c) fan, condition satisfied: {n} instances, deterministic != t0 in {bd}')
# (d)
mn=1
for _ in range(60):
    q=rng.uniform(.18,.19,5); b=rng.uniform(.039,.041,10)
    mn=min(mn,rstar(q,b)[0]-rand_val(q,b))
print(f'(d) open box: min (deterministic - randomized) over 60 random points = {mn:.4f} (claimed > 0.0145)   [{time.time()-T0:.0f}s]')
