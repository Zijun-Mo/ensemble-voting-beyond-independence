"""Independent checks for the B_graph NP-hardness report.
(1) J4° (9 vertices): B_graph via full judicious-partition MILP over all 2^8 complementary pairs (no gadget assumptions).
(2) Voting instance on K=9 in R0: q_i=1/9, b=1/36 on the 10 paid edges, 0 on the 4 vertical edges, cap (=1/9) elsewhere.
    Exact R*_det by the full no-big-M certificate MILP; compare with gamma*14/5 and with min(q_min, c*).
(3) Fractional-chromatic threshold theorem on further graphs (Petersen, Groetzsch, C7-complement, W5, K4) via an
    independently written reduced MILP (vertical edges free) and an LP for chi_f."""
import itertools, time, numpy as np, networkx as nx
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
T0=time.time()
def bgraph(nv, edges, cost):
    m=len(edges); reps=[S for S in range(1<<nv) if not S&1]; nz=len(reps); N=m+nz; R=[];lo=[]
    for t,S in enumerate(reps):
        A=np.zeros(N); D=np.zeros(N)
        for k,(u,v) in enumerate(edges):
            iu,iv=(S>>u)&1,(S>>v)&1
            if iu and iv: A[k]=1
            if not iu and not iv: D[k]=1
        r=A.copy(); r[m+t]=-1; R.append(r); lo.append(0)        # A_S >= z
        r=D.copy(); r[m+t]=1; R.append(r); lo.append(1)         # D_S >= 1-z
    res=milp(np.r_[cost,np.zeros(nz)],constraints=LinearConstraint(np.array(R),lo,np.inf),integrality=np.r_[np.zeros(m),np.ones(nz)],bounds=Bounds(0,1),options={'mip_rel_gap':1e-10})
    return res.fun
# J4°: vertices 0=B0, 1..4=T1..T4, 5..8=B1..B4
E=[(i,4+i) for i in range(1,5)]+[(0,4+i) for i in range(1,5)]+list(itertools.combinations(range(1,5),2))
cost=np.r_[np.zeros(4),np.ones(10)]
print(f'(1) J4° full-partition MILP: B = {bgraph(9,E,cost):.6f} (claim 14/5 = 2.8); c* = 3')
# (2) K=9 voting instance
K=9; Z=np.array([[(m>>i)&1 for i in range(K)] for m in range(1<<K)]); PAIRS=list(itertools.combinations(range(K),2))
q=np.full(K,1/9); b=np.full(len(PAIRS),1/9)
for (u,v),c in zip(E,cost): b[PAIRS.index(tuple(sorted((u,v))))]=c/36
F=np.hstack([np.ones((len(Z),1)),Z,np.array([Z[:,i]*Z[:,j] for i,j in PAIRS]).T]); nc=F.shape[1]
REPS=[m for m in range(1<<K) if not m&1]; nb=len(REPS); N=nc+nb; R=[];lo=[]
for t,m in enumerate(REPS):
    r=np.zeros(N); r[:nc]=F[m]; r[nc+t]=-1; R.append(r); lo.append(0)
    r=np.zeros(N); r[:nc]=F[(1<<K)-1-m]; r[nc+t]=1; R.append(r); lo.append(1)
lb=np.r_[-np.inf*np.ones(1+K),np.zeros(len(PAIRS)),np.zeros(nb)]; ub=np.r_[np.inf*np.ones(nc),np.ones(nb)]
res=milp(np.r_[1,q,b,np.zeros(nb)],constraints=LinearConstraint(np.array(R),lo,np.inf),integrality=np.r_[np.zeros(nc),np.ones(nb)],bounds=Bounds(lb,ub),options={'mip_rel_gap':1e-9,'time_limit':150})
print(f'(2) K=9 voting instance in R0 (sum q = 1, feasible via mutually exclusive errors): R*_det = {res.fun:.6f}; '
      f'gamma*14/5 = {14/180:.6f}; min(q_min, c*) = {min(1/9,3/36):.6f}   [{time.time()-T0:.0f}s]')
# (3)
def chif(G):
    n=G.number_of_nodes(); rows=[]
    for S in range(1,1<<n):
        vs=[v for v in range(n) if S>>v&1]
        if all(not G.has_edge(u,v) for u,v in itertools.combinations(vs,2)): rows.append([1 if v in vs else 0 for v in range(n)])
    return -linprog(-np.ones(n),A_ub=np.array(rows),b_ub=np.ones(len(rows)),bounds=(0,None),method='highs').fun
def BHG(G,beta):
    n=G.number_of_nodes(); Ed=list(G.edges); m=len(Ed); flex=[];fixed=[]
    for X in range(1<<n):
        (flex if any((X>>u)&1 and (X>>v)&1 for u,v in Ed) else fixed).append(X)
    N=n+m+len(flex); R=[];lo=[]
    def Arow(X):
        r=np.zeros(N)
        for v in range(n):
            if X>>v&1: r[v]=1
        for k,(u,v) in enumerate(Ed):
            if not (X>>u)&1 and not (X>>v)&1: r[n+k]=1
        return r
    for X in fixed: R.append(Arow(X)); lo.append(1)
    for t,X in enumerate(flex):
        r=Arow(X); r[n+m+t]=-1; R.append(r); lo.append(0)
        d=np.zeros(N)
        for k,(u,v) in enumerate(Ed):
            if (X>>u)&1 and (X>>v)&1: d[n+k]=1
        d[n+m+t]=1; R.append(d); lo.append(1)
    res=milp(np.r_[np.ones(n),beta*np.ones(m),np.zeros(len(flex))],constraints=LinearConstraint(np.array(R),lo,np.inf),integrality=np.r_[np.zeros(n+m),np.ones(len(flex))],bounds=Bounds(0,1),options={'mip_rel_gap':1e-10,'time_limit':40})
    return res.fun
graphs={'Petersen':nx.petersen_graph(),'Groetzsch':nx.mycielski_graph(4),'C7-complement':nx.complement(nx.cycle_graph(7)),'W5':nx.wheel_graph(6),'K4':nx.complete_graph(4)}
for name,G in graphs.items():
    G=nx.convert_node_labels_to_integers(G); cf=chif(G); out=[]
    for beta in (1.0,1.25,1.5,1.75,2.0):
        if time.time()-T0>280: break
        B=BHG(G,beta); out.append(f'beta={beta}: B={B:.4f}{"<" if B<2+beta-1e-6 else "="}{2+beta} {"OK" if (B<2+beta-1e-6)==(cf>2+beta+1e-9) else "MISMATCH"}')
    print(f'(3) {name} (n={G.number_of_nodes()}, chi_f={cf:.3f}): '+'; '.join(out))
print(f'[{time.time()-T0:.0f}s]')
