import numpy as np, itertools, time
exec(open('check_oddcycle.py').read().split('# (4)')[0].replace("print(f'(3)","pass #").replace("while time.time()-T0<120","while False"))
from scipy.optimize import milp, LinearConstraint, Bounds
G20=[];lb=[]
for i in range(5):
    r=np.zeros(20); r[[k for k,(u,v) in enumerate(EDGES) if i not in (u,v)]]=1; G20.append(r); lb.append(1)
for k,(i,j) in enumerate(EDGES):
    r=np.zeros(20); r[k]=1; r[10+k]=-1; G20.append(r); lb.append(0)
    r=np.zeros(20); r[10+k]=1; r[[m for m,(u,v) in enumerate(EDGES) if len({u,v,i,j})==4]]=1; G20.append(r); lb.append(1)
G20=np.array(G20)
def g20(b): return milp(np.r_[b,np.zeros(10)],integrality=np.r_[np.zeros(10),np.ones(10)],bounds=Bounds(0,1),constraints=LinearConstraint(G20,lb,np.inf),options={'mip_rel_gap':1e-11}).fun
rng=np.random.default_rng(3); T0=time.time(); shown=0; n=0; bad_pure=0; bad20=0
while time.time()-T0<150:
    b=rng.uniform(0,1,10)**rng.choice([1,3,6]); b*=(rng.uniform(size=10)<rng.choice([0.6,1.0]))
    c=cst(b); p=bgraph(b); t=g20(b); n+=1
    if abs(p-c)>1e-6: bad_pure+=1
    if abs(t-c)>1e-6: bad20+=1
    if (abs(p-c)>1e-6 or abs(t-c)>1e-6) and shown<3:
        shown+=1; print('b',np.round(b,4),'c*',c,'pureMILP',p,'20varMILP',t)
print(n,'instances; pure-MILP mismatches',bad_pure,'; 20-var MILP mismatches',bad20)
