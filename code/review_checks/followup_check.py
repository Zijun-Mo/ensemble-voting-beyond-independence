import numpy as np, itertools, time
exec(open('graph_cert_exactness.py').read().split("for K in [3,4,5,6,7]:")[0])
T0=time.time()
# (1) the K=5 symmetric counterexample and the positive-width family
K=5; Z=states(K); pairs=list(itertools.combinations(range(K),2))
for qv,bv in [(0.3,0.05),(0.3,0.051),(0.3,0.06),(0.3,0.0666),(0.25,0.03),(0.35,0.08)]:
    q=np.full(K,qv); b={p:bv for p in pairs}
    full=solve(K,Z,q,b,pairs,False); Bg=solve(K,Z,q,b,pairs,True)
    print(f"q={qv}, b={bv}: R*_det={full:.4f}  B_graph={Bg:.4f}  min(q,Bg)={min(qv,Bg):.4f}  moment 1-5q+10b={1-5*qv+10*bv:.4f}  in R0(sum q<=1)? {5*qv<=1}")
# (2) extended candidate: min{q_min, B_graph, shifted-moment certificates on odd subsets T}
#     majority on T (|T|=2r+1) fails => N_T>=r+1;  1{N>=r+1} <= (N-t)(N-t+1)/2 for integer t<=r  (t=r is tightest-looking; scan t=1..r)
def moment_bound(q,b,T):
    r=(len(T)-1)//2; best=1
    sq=sum(q[i] for i in T); sb=sum(b[(i,j)] for i,j in itertools.combinations(sorted(T),2))
    for t in range(1,r+1):  # E[(N-t)(N-t+1)/2] = E C(N,2) - (t-1)E N + t(t-1)/2
        best=min(best, sb-(t-1)*sq+t*(t-1)/2)
    return max(best,0)
rng=np.random.default_rng(21)
for K in [5,6,7]:
    Z=states(K); pairs=list(itertools.combinations(range(K),2)); n=g1=g2=0; gaps=[]
    while time.time()-T0<70*(K-4) and n<150:
        q=rng.uniform(0.2,0.45,K)
        if not (q.sum()>1 and q.sum()-(K-2)*q.min()<=1): continue   # sample R \ R0 only
        b={p:min(q[p[0]],q[p[1]])*rng.uniform(0,0.6) for p in pairs}
        full=solve(K,Z,q,b,pairs,False)
        if not np.isfinite(full): continue
        n+=1; Bg=solve(K,Z,q,b,pairs,True); c1=min(q.min(),Bg)
        mom=min(moment_bound(q,b,T) for s in range(3,K+1,2) for T in itertools.combinations(range(K),s))
        c2=min(c1,mom)
        g1+= c1>full+1e-9; g2+= c2>full+1e-9
        if c2>full+1e-9: gaps.append(c2-full)
    print(f"K={K}, R\\R0 samples={n}: graph-cert fails {g1}; graph+moment fails {g2}; max remaining gap={max(gaps) if gaps else 0:.4f}; t={time.time()-T0:.0f}s",flush=True)
