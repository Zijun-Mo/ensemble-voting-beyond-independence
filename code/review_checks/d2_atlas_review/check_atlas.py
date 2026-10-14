"""Independent checks of the atlas report: (1) 20-var graph MILP vs pure-quadratic certificate MILP;
(2) zero-intercept theorem; (3) atlas formula and positive-intercept corollary vs ground truth
(min over all 1882 signed threshold rules of the fixed-rule worst-case LP)."""
import numpy as np, itertools, json, time
from scipy.optimize import milp, LinearConstraint, Bounds
from exact_core import Z, FEAT, EDGES, K, P, REPS, fixed_risk, feasible
T0=time.time(); rng=np.random.default_rng(7)
# --- report's 20-variable graph MILP
G=[];lb=[]
for i in range(5):
    r=np.zeros(20); r[[k for k,(u,v) in enumerate(EDGES) if i not in (u,v)]]=1; G.append(r); lb.append(1)
for k,(i,j) in enumerate(EDGES):
    r=np.zeros(20); r[k]=1; r[10+k]=-1; G.append(r); lb.append(0)
    r=np.zeros(20); r[10+k]=1; r[[m for m,(u,v) in enumerate(EDGES) if len({u,v,i,j})==4]]=1; G.append(r); lb.append(1)
G=np.array(G)
def graph20(b):
    r=milp(np.r_[b,np.zeros(10)],integrality=np.r_[np.zeros(10),np.ones(10)],bounds=Bounds(0,1),
           constraints=LinearConstraint(G,lb,np.inf),options={'mip_rel_gap':1e-11}); return r.fun
# --- certificate MILP with optional restrictions
def certmilp(q,b,mode):
    nc=16; nb=16; n=nc+nb; rows=[];lo=[]
    for t,m in enumerate(REPS):
        a=np.zeros(n); a[:nc]=FEAT[m]; a[nc+t]=-1; rows.append(a); lo.append(0)
        a=np.zeros(n); a[:nc]=FEAT[31-m]; a[nc+t]=1; rows.append(a); lo.append(1)
    lbv=np.r_[-np.inf*np.ones(6),np.zeros(10),np.zeros(nb)]; ubv=np.r_[np.inf*np.ones(nc),np.ones(nb)]
    if mode in ('a0','pure'): lbv[0]=ubv[0]=0
    if mode=='pure': lbv[1:6]=ubv[1:6]=0
    r=milp(np.r_[1,q,b,np.zeros(nb)],constraints=LinearConstraint(np.array(rows),lo,np.inf),
           integrality=np.r_[np.zeros(nc),np.ones(nb)],bounds=Bounds(lbv,ubv),options={'mip_rel_gap':1e-11}); return r.fun
def sample(region):
    while True:
        k=rng.integers(2,12); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k))
        q=Z.T@pi
        if q.min()<0.02 or q.max()>=0.5: continue
        inR=q.sum()-3*q.min()<=1
        if (region=='R\\R0' and not (inR and q.sum()>1)) or (region=='outR' and inR): continue
        J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
        return q, J+(cap-J)*rng.uniform(0,1,P)**2*(rng.uniform(size=P)<0.5)
# (1)
d1=0
for _ in range(300):
    b=rng.uniform(0,0.3,10)*(rng.uniform(size=10)<0.8); d1=max(d1,abs(graph20(b)-certmilp(np.full(5,.3),b,'pure')))
print(f'(1) graph20 vs pure-quadratic MILP, 300 random b: max |diff| = {d1:.2e}   [{time.time()-T0:.0f}s]')
# (2)
d2=0; cnt=0
for reg in ('R\\R0','outR'):
    for _ in range(100):
        q,b=sample(reg); d2=max(d2,abs(certmilp(q,b,'a0')-min(q.min(),graph20(b)))); cnt+=1
print(f'(2) zero-intercept optimum vs min(q_min,B_graph), {cnt} instances (R\\R0 and outside R): max |diff| = {d2:.2e}   [{time.time()-T0:.0f}s]')
# (3) atlas
data=json.load(open('rule_vertices.json'))
PERMS=list(itertools.permutations(range(5)))
def idx(p):
    ix=[0]+[0]*15
    for i in range(5): ix[1+p[i]]=1+i
    for e,(i,j) in enumerate(EDGES): ix[6+EDGES.index(tuple(sorted((p[i],p[j]))))]=6+e
    return ix
V0=np.array([[float(eval(x)) for x in v] for d in data for v in d['vertices']])
ALL=np.vstack([V0[:,idx(p)] for p in PERMS]); POS=ALL[ALL[:,0]>1e-12]
print(f'    atlas: {len(V0)} orbit-rep vertices -> {len(ALL)} with images; a>0 part {len(POS)}')
rules=set()
for w0 in [(1,0,0,0,0),(3,1,1,1,1),(2,1,1,1,0),(3,2,2,1,1),(1,1,1,0,0),(2,2,1,1,1),(1,1,1,1,1)]:
    for w in set(itertools.permutations(w0)):
        for sg in itertools.product((1,-1),repeat=5):
            rules.add(tuple(int(sum(w[i]*sg[i]*(2*z[i]-1) for i in range(5))>0) for z in Z))
rules=[np.array(f,float) for f in rules]; assert len(rules)==1882
d3=d4=0; n3=0
for reg in ('R\\R0','outR'):
    for _ in range(9):
        if time.time()-T0>270: break
        q,b=sample(reg); th=np.r_[1,q,b]
        truth=min(fixed_risk(f,q,b) for f in rules)
        atlas=(ALL@th).min(); cor=min(q.min(),graph20(b),(POS@th).min())
        d3=max(d3,abs(truth-atlas)); d4=max(d4,abs(truth-cor)); n3+=1
print(f'(3) {n3} instances: max|truth-atlas| = {d3:.2e}, max|truth-corollary| = {d4:.2e}   [{time.time()-T0:.0f}s]')
