import numpy as np, time, itertools
exec(open('graph_cert_exactness.py').read().split("for K in [3,4,5,6,7]:")[0])
T0=time.time(); rng=np.random.default_rng(9)
for K in [3,4,5,6,7]:
    Z=states(K); pairs=list(itertools.combinations(range(K),2)); n=bad=nout=0; gaps=[]; ex=None
    while time.time()-T0<50*(K-2) and n<250:
        q=rng.uniform(0.15,0.49,K)
        b={p:min(q[p[0]],q[p[1]])*(rng.uniform(0,0.9) if rng.random()<0.8 else 1.0) for p in pairs}
        full=solve(K,Z,q,b,pairs,False)
        if not np.isfinite(full): continue
        cand=min(q.min(),solve(K,Z,q,b,pairs,True)); n+=1; nout+= q.sum()-(K-2)*q.min()>1
        if cand>full+1e-9:
            bad+=1; gaps.append(cand-full)
            if ex is None: ex=(np.round(q,3),{k:round(v,3) for k,v in b.items()},round(full,4),round(cand,4))
    print(f"K={K}: n={n} (outside R: {nout}); min{{q_min,B_graph}} > R* in {bad}; max gap={max(gaps) if gaps else 0:.4f}; t={time.time()-T0:.0f}s",flush=True)
    if ex: print("   example:",ex)
