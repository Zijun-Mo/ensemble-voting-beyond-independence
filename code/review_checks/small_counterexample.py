"""Search explicit small counterexamples to T2/T3 on small triangle-free graphs H.
q_i = 1/(2K); b_e = beta*q on edges of H, b = q (non-binding) elsewhere.  c* = min(q, oddgirth*beta*q).
Exact worst-case risk of a given weighted vote by LP over 2^K error patterns (only edge constraints bind).
Rules tried: degree-weighted vote and simple majority, generic positive tie-breaks (several seeds)."""
import itertools, time, numpy as np, networkx as nx
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack
T0=time.time(); rng=np.random.default_rng(11)
def graphs():
    yield "C5", nx.cycle_graph(5)
    yield "Wagner(8)", nx.circulant_graph(8,[1,4])
    yield "Petersen(10)", nx.petersen_graph()
    yield "Grotzsch(11)", nx.mycielski_graph(4)
    yield "Andrasfai4(11)", nx.circulant_graph(11,[1,4,7,10])
    yield "Clebsch(16)", nx.circulant_graph(16,[1,4,6,7]) if False else nx.from_edgelist([(a,b) for a in range(16) for b in range(a+1,16) if bin(a^b).count('1') in (1,4)])
def odd_girth(G):
    best=99
    for s in G.nodes:
        dist={s:0}; par={s:None}; Q=[s]
        for u in Q:
            for v in G[u]:
                if v not in dist: dist[v]=dist[u]+1; Q.append(v)
                elif dist[v]==dist[u]: best=min(best,2*dist[u]+1)
    return best
def worst_risk(f,Z,q,edges,b):
    n=len(Z)
    A_eq=vstack([csr_matrix(np.ones((1,n)))]+[csr_matrix(Z[:,i][None,:]) for i in range(Z.shape[1])])
    A_ub=vstack([csr_matrix((Z[:,i]*Z[:,j])[None,:]) for i,j in edges])
    r=linprog(-f,A_ub=A_ub,b_ub=np.full(len(edges),b),A_eq=A_eq,b_eq=np.r_[1,q],bounds=(0,None),method='highs')
    return -r.fun
for name,G in graphs():
    if time.time()-T0>230: print("time budget reached"); break
    G=nx.convert_node_labels_to_integers(G); K=G.number_of_nodes(); E=list(G.edges())
    assert sum(nx.triangles(G).values())==0
    og=odd_girth(G); deg=np.array([G.degree(i) for i in range(K)],float)
    Z=np.array(list(itertools.product([0,1],repeat=K)),dtype=np.int8); S=(2*Z-1).astype(float)
    q=np.full(K,1/(2*K)); out=[]
    for beta in [0.08,0.12,0.16,0.20,0.25]:
        b=beta*q[0]; cstar=min(q[0],og*b); best=9
        for base in (deg,np.ones(K)):
            for s in range(3):
                w=base+rng.uniform(1e-6,1e-5,K); f=(S@w>0).astype(float)
                best=min(best,worst_risk(f,Z,q,E,b))
        out.append((beta, round(best/q[0],4), round(cstar/q[0],4), best<cstar-1e-12))
    print(f"{name}: K={K}, edges={len(E)}, odd girth={og}, maxcut/m≈{nx.algorithms.approximation.maxcut.one_exchange(G,seed=1)[0]/len(E):.3f} (local)")
    for beta,r,c,flag in out: print(f"   beta={beta}: best risk/q={r}  c*/q={c}  counterexample={flag}")
    print(f"   t={time.time()-T0:.0f}s",flush=True)
