"""Independent checks of the 808-graph / chorded-cycle report.
(a) K=4 mixed (randomized) dual {Q>=0, Q(x)+Q(1-x)>=1, gamma>=0}: exact vertex count; randomized value == 12-term == deterministic
(b) K=5: random parameters whose effective graph lies in the 808 solved classes: exact R* vs graph formula
(c) C5+chord: random parameters with the triangle condition: exact R* vs t; and how often t is still exact without it
(d) the report's boundary point (q=0.45): exact R*"""
import itertools, time, numpy as np, networkx as nx
from fractions import Fraction as Fr
import cdd, cdd.gmp as G
from scipy.optimize import linprog
from exact_core import Z, EDGES, rstar, fixed_risk
T0=time.time(); rng=np.random.default_rng(9)
# (a)
E4=list(itertools.combinations(range(4),2)); X4=[[(m>>i)&1 for i in range(4)] for m in range(16)]
F4=[[1,*x,*(x[i]*x[j] for i,j in E4)] for x in X4]
rows=[[0]+F4[m] for m in range(16)]+[[-1]+[a+b for a,b in zip(F4[m],F4[15-m])] for m in range(16) if m<15-m]+[[0]*6+[int(k==e) for k in range(6)] for e in range(6)]
g=G.copy_generators(G.polyhedron_from_matrix(G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY)))
VM=np.array([[float(Fr(x)) for x in r[1:]] for r in g.array if r[0]==1])
C12=[]
for i in range(4): c=np.zeros(11); c[1+i]=1; C12.append(c)
for T in itertools.combinations(range(4),3):
    c=np.zeros(11)
    for e in itertools.combinations(T,2): c[5+E4.index(e)]=1
    C12.append(c)
for r in range(4): C12.append(np.r_[1,-np.ones(4),[2 if r in e else 1 for e in E4]])
C12=np.array(C12); Z4=np.array(X4); P4=np.array([Z4[:,i]*Z4[:,j] for i,j in E4])
mx=0
for _ in range(400):
    k=rng.integers(2,10); S=rng.choice(16,k,replace=False); pi=np.zeros(16); pi[S]=rng.dirichlet(np.ones(k)); q=Z4.T@pi
    if q.max()>=.5: continue
    J=P4@pi; cap=np.array([min(q[i],q[j]) for i,j in E4]); b=J+(cap-J)*rng.uniform(0,1,6)**2*(rng.uniform(size=6)<.5)
    th=np.r_[1,q,b]; mx=max(mx,abs((VM@th).min()-(C12@th).min()))
print(f'(a) K=4 mixed dual: {len(VM)} vertices (report: 24); randomized value vs 12-term, max diff {mx:.2e}')
# helpers
def odd_cycles(Hg):
    out=[]
    for k in (3,5):
        for S in itertools.combinations(range(5),k):
            for perm in itertools.permutations(S[1:]):
                if k==5 and perm[0]>perm[-1]: continue
                if k==3 and perm[0]>perm[1]: continue
                cyc=(S[0],)+perm; es=[tuple(sorted((cyc[i],cyc[(i+1)%k]))) for i in range(k)]
                if all(Hg.has_edge(*e) for e in es): out.append(es)
    return out
def formula(q,b,Hg):
    bd={e:b[EDGES.index(e)] for e in EDGES}; terms=[q.min()]
    terms+=[sum(bd[e] for e in C) for C in odd_cycles(Hg)]
    for B in itertools.combinations(range(5),4):
        if all(Hg.has_edge(*e) for e in itertools.combinations(B,2)):
            for r in B: terms.append(1-sum(q[i] for i in B)+sum(bd[e] for e in itertools.combinations(B,2) if r not in e)+2*sum(bd[tuple(sorted((r,j)))] for j in B if j!=r))
    return min(terms)
def category(Hg):
    if Hg.number_of_edges()==5 and all(d==2 for _,d in Hg.degree) and nx.is_connected(Hg): return 'C5'
    if not any(len(B)==5 for B in nx.biconnected_components(Hg)): return 'blocks<=4'
    return 'bipartite5' if nx.is_bipartite(Hg) else 'remaining'
RULES=None
def truth_rstar(q,b):
    global RULES
    if RULES is None:
        R=set()
        for w0 in [(1,0,0,0,0),(3,1,1,1,1),(2,1,1,1,0),(3,2,2,1,1),(1,1,1,0,0),(2,2,1,1,1),(1,1,1,1,1)]:
            for w in set(itertools.permutations(w0)):
                for sg in itertools.product((1,-1),repeat=5): R.add(tuple(int(sum(w[i]*sg[i]*(2*z[i]-1) for i in range(5))>0) for z in Z))
        RULES=[np.array(f,float) for f in R]
    return min(fixed_risk(f,q,b) for f in RULES)
def params_for(H,force=None):
    while True:
        k=rng.integers(2,12); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k)*rng.choice([.5,1,3])); q=Z.T@pi
        if q.min()<0.02 or q.max()>=.5: continue
        J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
        b=np.where([H.has_edge(*e) for e in EDGES], J+(cap-J)*rng.uniform(0,.9,10)**2, cap)
        Heff=nx.Graph(); Heff.add_nodes_from(range(5)); Heff.add_edges_from(e for k2,e in enumerate(EDGES) if b[k2]<cap[k2]-1e-12)
        return q,b,Heff
# (b)
masks={'blocks<=4':[],'C5':[],'bipartite5':[]}
for mask in range(1024):
    Hg=nx.Graph(); Hg.add_nodes_from(range(5)); Hg.add_edges_from(e for j,e in enumerate(EDGES) if mask>>j&1)
    c=category(Hg)
    if c in masks: masks[c].append(Hg)
stat={}; bad=[]
while time.time()-T0<170:
    c=rng.choice(list(masks)); H=masks[c][rng.integers(len(masks[c]))]
    q,b,He=params_for(H); ce=category(He)
    if ce=='remaining': continue
    r=rstar(q,b)[0]; f=formula(q,b,He)
    s=stat.setdefault(ce,[0,0,0]); s[0]+=1; s[2]+=q.sum()>1
    if abs(r-f)>1e-7:
        t=truth_rstar(q,b)
        if abs(t-f)>1e-7: s[1]+=1; bad.append((ce,q,b,t,f))
print('(b) solved classes: '+'; '.join(f'{k}: {v[0]} instances ({v[2]} with sum q>1), mismatches {v[1]}' for k,v in stat.items()))
# (c)
Hc=nx.Graph(); Hc.add_nodes_from(range(5)); Hc.add_edges_from([(0,1),(1,2),(2,3),(3,4),(0,4),(0,2)])
nin=bin_=nout=bout=0
while time.time()-T0<270:
    q,b,He=params_for(Hc)
    if set(map(tuple,map(sorted,He.edges)))!={(0,1),(1,2),(2,3),(3,4),(0,4),(0,2)}: continue
    bd={e:b[EDGES.index(e)] for e in EDGES}
    t=min(q.min(),bd[(0,1)]+bd[(1,2)]+bd[(0,2)],bd[(0,1)]+bd[(1,2)]+bd[(2,3)]+bd[(3,4)]+bd[(0,4)])
    r=rstar(q,b)[0]
    if t>=q[0]+q[1]+q[2]-1:
        nin+=1
        if abs(r-t)>1e-7 and abs(truth_rstar(q,b)-t)>1e-7: bin_+=1
    else:
        nout+=1; bout+= abs(r-t)>1e-7
print(f'(c) C5+chord: {nin} instances with the condition, mismatches {bin_}; without condition: {nout} instances, t != R* in {bout}')
# (d)
q=np.full(5,.45); b=np.array([.05 if e in [(0,1),(1,2),(2,3),(3,4),(0,4)] else (.40 if e==(0,2) else .45) for e in EDGES])
print(f'(d) boundary point q=0.45: t=0.25, exact R* (MILP) = {rstar(q,b)[0]:.6f}, by 1882-rule LPs = {truth_rstar(q,b):.6f}   [{time.time()-T0:.0f}s]')
