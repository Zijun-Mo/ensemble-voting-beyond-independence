"""Task-1 report checks (K=5): (0) exact Lipschitz constants from the 17,508 dual vertices;
(1) coverage + gap of the graph certificate C_G = min{q+_min, c*(b_hat)} (K=5: B_graph=c*);
(2) full-ambiguity certificate C_full via 1,882 threshold-rule LPs, inside and outside R0."""
import numpy as np, json, itertools, time
from scipy.optimize import linprog
from exact_core import Z, EDGES, rstar
T0=time.time(); rng=np.random.default_rng(11); K=5; m=10
V=np.array([[float(eval(x)) for x in v] for d in json.load(open('rule_vertices.json')) for v in d['vertices']])
L5=np.max(np.abs(V[:,1:6]).sum(1)+V[:,6:].sum(1))
print(f'(0) exact L_5 over all dual vertices = {L5:g}; crude bound (K+m)d^(d/2) = {(K+m)*16**8:.2e}; zero-intercept part: C_0 = min of 27 atoms, slope <= 5')
tri=[tuple(EDGES.index(e) for e in itertools.combinations(S,2)) for S in itertools.combinations(range(5),3)]
pent=[]
for r in itertools.permutations(range(1,5)):
    if r[0]>r[-1]: continue
    c=(0,)+r; pent.append(tuple(sorted(EDGES.index(tuple(sorted((c[i],c[(i+1)%5])))) for i in range(5))))
CYC=tri+pent
rules=set()
for w0 in [(1,0,0,0,0),(3,1,1,1,1),(2,1,1,1,0),(3,2,2,1,1),(1,1,1,0,0),(2,2,1,1,1),(1,1,1,1,1)]:
    for w in set(itertools.permutations(w0)):
        for sg in itertools.product((1,-1),repeat=5):
            rules.add(tuple(int(sum(w[i]*sg[i]*(2*z[i]-1) for i in range(5))>0) for z in Z))
RULES=[np.array(f,float) for f in rules]
PAIR=np.array([Z[:,i]*Z[:,j] for i,j in EDGES])
def estimates(pi0,n,delta):
    H=np.log(2*(K+m)/delta); eps=np.sqrt(H/(2*n))
    lab=Z[rng.choice(32,n,p=pi0)]; unl=Z[rng.choice(32,n,p=pi0)]
    qh=lab.mean(0); Dh=np.array([(unl[:,i]!=unl[:,j]).mean() for i,j in EDGES])
    qp=np.clip(qh+eps,0,.5); qm=np.clip(qh-eps,0,.5); Dm=np.clip(Dh-eps,0,1)
    bh=np.array([min(qp[i],qp[j],np.clip((qp[i]+qp[j]-Dm[k])/2,0,1)) for k,(i,j) in enumerate(EDGES)])
    return qp,qm,bh,eps
def graph_cert(qp,bh):
    cyc=[sum(bh[e] for e in C) for C in CYC]; k=int(np.argmin(cyc)); i=int(np.argmin(qp))
    if qp[i]<=cyc[k]: f=Z[:,i].astype(float); return qp[i],f
    verts=sorted({v for e in CYC[k] for v in EDGES[e]})
    f=(Z[:,verts].sum(1)*2>len(verts)).astype(float); return cyc[k],f
def full_cert(qp,qm,bh):
    Aub=np.vstack([Z.T,-Z.T,PAIR]); bub=np.r_[qp,-qm,bh]
    best=1.0
    for f in RULES:
        r=linprog(-f,A_ub=Aub,b_ub=bub,A_eq=np.ones((1,32)),b_eq=[1],bounds=(0,None),method='highs')
        if r.status==0: best=min(best,-r.fun)
    return best
# population laws
piA=np.zeros(32)   # inside R0: sparse random law
while True:
    S=rng.choice(32,8,replace=False); piA[:]=0; piA[S]=rng.dirichlet(np.ones(8)); q=Z.T@piA
    if q.sum()<=1 and q.max()<.5 and q.min()>.05: break
piB=np.zeros(32)   # outside R0: q=0.3, e=0.05 (all 2-error patterns 0.05, single errors 0.1)
for mm in range(32):
    c=bin(mm).count('1'); piB[mm]=0.05 if c==2 else (0.1 if c==1 else 0)
for name,pi0 in (('R0 law',piA),('q=.3,e=.05 law',piB)):
    q=Z.T@pi0; e=PAIR@pi0; V0=rstar(q,e)[0]
    print(f'\n[{name}] sum q={q.sum():.3f}, V0={V0:.4f}, C_0(q,e)=min(q_min,c*)={min(q.min(),min(sum(e[x] for x in C) for C in CYC)):.4f}')
    for n in (500,5000,50000):
        cov=0; gaps=[]; bnd=None
        for rep in range(200):
            qp,qm,bh,eps=estimates(pi0,n,0.1); c,f=graph_cert(qp,bh)
            cov+= (f@pi0<=c+1e-12); gaps.append(c-V0); bnd=5*(3*eps)  # L0=5, 2eps_l+eps_u=3eps
        print(f'  n={n:6d}: graph cert coverage {cov}/200, gap C_G-V0 mean {np.mean(gaps):.4f} max {np.max(gaps):.4f} (Thm A bound in R0: {bnd:.4f})')
    for n in (2000,200000):
        if time.time()-T0>240: break
        qp,qm,bh,eps=estimates(pi0,n,0.1); cf=full_cert(qp,qm,bh)
        print(f'  n={n:6d}: C_full={cf:.4f} (gap {cf-V0:.4f}; Thm B bound {L5*3*eps:.3f})')
print(f'[{time.time()-T0:.0f}s]')
