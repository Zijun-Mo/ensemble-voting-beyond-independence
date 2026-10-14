import pickle, itertools, numpy as np, time
from scipy.optimize import milp, LinearConstraint, Bounds
u7,uA,uB=pickle.load(open('c5k2.pkl','rb'))
def bgraph(n,u):
    E=list(itertools.combinations(range(n),2)); cost=np.array([float(u[e]) for e in E])
    m=len(E); reps=[S for S in range(1<<n) if not S&1]; nz=len(reps); N=m+nz; R=[];lo=[]
    for t,S in enumerate(reps):
        a=np.zeros(N); d=np.zeros(N)
        for k,(x,y) in enumerate(E):
            if (S>>x)&1 and (S>>y)&1: a[k]=1
            if not (S>>x)&1 and not (S>>y)&1: d[k]=1
        a[m+t]=-1; R.append(a); lo.append(0); d[m+t]=1; R.append(d); lo.append(1)
    r=milp(np.r_[cost,np.zeros(nz)],constraints=LinearConstraint(np.array(R),lo,np.inf),integrality=np.r_[np.zeros(m),np.ones(nz)],bounds=Bounds(0,1),options={'mip_rel_gap':1e-12,'time_limit':100})
    return r.fun, getattr(r,'mip_dual_bound',None), r.status
for name,n,u in (('C5+K2',7,u7),('LiftA',8,uA),('LiftB',8,uB)):
    t=time.time(); print(name, bgraph(n,u), f'{time.time()-t:.1f}s')
