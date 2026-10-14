import numpy as np, itertools, time
exec(open('graph_cert_exactness.py').read().split("for K in [3,4,5,6,7]:")[0])
T0=time.time(); rng=np.random.default_rng(33)
def cstar(K,b):
    best=np.inf
    for L in range(3,K+1,2):
        for cyc in itertools.permutations(range(K),L):
            if cyc[0]!=min(cyc) or cyc[1]>cyc[-1]: continue
            best=min(best,sum(b[tuple(sorted((cyc[k],cyc[(k+1)%L])))] for k in range(L)))
    return best
for K in [5,6,7]:
    Z=states(K); pairs=list(itertools.combinations(range(K),2)); n=strict=0; best_gap=0; ex=None
    while time.time()-T0<80*(K-4) and n<300:
        q=rng.uniform(0.02,0.25,K); q=q/q.sum()*rng.uniform(0.5,1.0)        # sum q <= 1 (R0)
        mode=rng.random()
        if mode<0.5:   # structured: low co-error on a sparse random graph, high elsewhere
            b={p:min(q[p[0]],q[p[1]])*(rng.uniform(0,0.25) if rng.random()<0.5 else 1.0) for p in pairs}
        else:
            b={p:min(q[p[0]],q[p[1]])*rng.uniform(0,1) for p in pairs}
        Bg=solve(K,Z,q,b,pairs,True); c=cstar(K,b); n+=1
        if Bg<min(q.min(),c)-1e-9:
            strict+=1; g=min(q.min(),c)-Bg
            if g>best_gap: best_gap=g; ex=(np.round(q,4),{k:round(v,4) for k,v in b.items()},round(Bg,4),round(c,4))
    print(f"K={K}: R0 samples={n}; B_graph < min(q_min,c*) in {strict}; largest gap={best_gap:.4f}; t={time.time()-T0:.0f}s",flush=True)
    if ex: print("   example:",ex)
