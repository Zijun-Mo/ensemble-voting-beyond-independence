"""Independent checks of the K=4 closed form / isolated-core report (all runs < 5 min).
(a) exact (GMP) vertex enumeration of the 14 K=4 signed-threshold-rule orbits; compare with the report's table
(b) every K=4 dual vertex is dominated by the 12 terms on the whole feasible domain (float LP, max Delta)
(c) R*_4 (min over 104 rules of the worst-case LP) == 12-term formula on random feasible instances
(d) K=5 isolated-vertex lemma: MILP with gamma_{0j}=0 forced  ==  min{q_0, R*_4(rest)}
(e) K=5 bipartite-support certificates never beat q_min (MILP with gamma restricted to a cut)
(f) graph counts 376 / 10 / 95 / 543 and 17 unlabeled spanning non-bipartite graphs
(g) anchor witness p=2/5: exact feasibility, H=0, all odd cycles > 0, MILP R*_5 = 0"""
import itertools, time, numpy as np
from fractions import Fraction as Fr
import cdd, cdd.gmp as G
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
T0=time.time(); rng=np.random.default_rng(5)
def setup(K):
    E=list(itertools.combinations(range(K),2)); Z=np.array([[(m>>i)&1 for i in range(K)] for m in range(2**K)])
    FE=np.array([np.r_[1,z,[z[i]*z[j] for i,j in E]] for z in Z],float); return E,Z,FE
E4,Z4,F4=setup(4); E5,Z5,F5=setup(5)
def rulef(w,Z): return tuple(int(sum(w[i]*(2*z[i]-1) for i in range(len(w)))>0) for z in Z)
# rules K=4
R4=set()
for w0 in [(1,0,0,0),(1,1,1,0),(2,1,1,1)]:
    for w in set(itertools.permutations(w0)):
        for sg in itertools.product((1,-1),repeat=4): R4.add(rulef(tuple(w[i]*sg[i] for i in range(4)),Z4))
assert len(R4)==104
def perm_rule(f,p,K):
    out=[0]*(2**K)
    for m in range(2**K): out[sum(((m>>i)&1)<<p[i] for i in range(K))]=f[m]
    return tuple(out)
orbs=[];seen=set()
for f in sorted(R4):
    if f in seen: continue
    o={perm_rule(f,p,4) for p in itertools.permutations(range(4))}; seen|=o; orbs.append((f,len(o)))
def vertices(f,FE,P):
    rows=[[-f[m]]+[int(v) for v in FE[m]] for m in range(len(f))]+[[0]*(1+1+FE.shape[1]-1-P)+[1 if k==e else 0 for k in range(P)] for e in range(P)]
    g=G.copy_generators(G.polyhedron_from_matrix(G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY)))
    return [tuple(Fr(x) for x in r[1:]) for r in g.array if r[0]==1]
VERTS={}
for f,sz in orbs: VERTS[f]=vertices(f,F4,6)
cnt=sorted(len(v) for v in VERTS.values())
print(f'(a) K=4: {len(orbs)} orbits; vertex counts per orbit {cnt}; total {sum(cnt)} (report: 262)')
# 12 terms
C12=[]
for i in range(4): c=np.zeros(11); c[1+i]=1; C12.append(c)
for T in itertools.combinations(range(4),3):
    c=np.zeros(11)
    for e in itertools.combinations(T,2): c[5+E4.index(e)]=1
    C12.append(c)
for r in range(4): C12.append(np.r_[1,-np.ones(4),[2 if r in e else 1 for e in E4]])
C12=np.array(C12)
def rowsA(cs,Z,P):
    A=np.zeros((len(cs),len(Z)+P))
    for k,c in enumerate(cs): c=np.array([float(x) for x in c]); A[k,:len(Z)]=c[0]+Z@c[1:1+Z.shape[1]]; A[k,len(Z):]=c[1+Z.shape[1]:]
    return A
AC=rowsA(C12,Z4,6); n=16+6+1; Aub=[];bub=[]
for e,(i,j) in enumerate(E4):
    r=np.zeros(n); r[:16]=Z4[:,i]*Z4[:,j]; r[16+e]=-1; Aub.append(r); bub.append(0)
    for k in (i,j): r=np.zeros(n); r[:16]=-Z4[:,k]; r[16+e]=1; Aub.append(r); bub.append(0)
for k in range(4): r=np.zeros(n); r[:16]=Z4[:,k]; Aub.append(r); bub.append(.5)
Aub=np.array(Aub); bub=np.array(bub); c=np.zeros(n); c[-1]=-1
mx=-1
for f,V in VERTS.items():
    for v in V:
        Av=rowsA([v],Z4,6)[0]; M=np.zeros((12,n)); M[:,:22]=-(AC-Av); M[:,-1]=1
        r=linprog(c,A_ub=np.vstack([Aub,M]),b_ub=np.r_[bub,np.zeros(12)],A_eq=np.r_[np.ones(16),np.zeros(7)][None],b_eq=[1],bounds=[(0,None)]*22+[(None,1)],method='highs')
        mx=max(mx,-r.fun)
print(f'(b) max Delta over all {sum(cnt)} orbit-rep vertices = {mx:.2e} (<=0 means dominated)')
# (c)
PAIR4=np.array([Z4[:,i]*Z4[:,j] for i,j in E4]); RL4=[np.array(f,float) for f in R4]
def rstar4(q,b):
    best=1
    for f in RL4:
        r=linprog(-f,A_ub=PAIR4,b_ub=b,A_eq=np.vstack([np.ones(16),Z4.T]),b_eq=np.r_[1,q],bounds=(0,None),method='highs'); best=min(best,-r.fun)
    return best
def form4(q,b): return (C12@np.r_[1,q,b]).min()
def sample(K,Z,E,cap_q=0.5):
    while True:
        k=rng.integers(2,10); S=rng.choice(len(Z),k,replace=False); pi=np.zeros(len(Z)); pi[S]=rng.dirichlet(np.ones(k)); q=Z.T@pi
        if q.min()>0.01 and q.max()<cap_q: break
    J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in E]); cap=np.array([min(q[i],q[j]) for i,j in E])
    return q, J+(cap-J)*rng.uniform(0,1,len(E))**2*(rng.uniform(size=len(E))<0.5)
d=0; nn=0; nout=0
while time.time()-T0<110:
    q,b=sample(4,Z4,E4); d=max(d,abs(rstar4(q,b)-form4(q,b))); nn+=1; nout+=q.sum()>1
print(f'(c) K=4: {nn} instances ({nout} with sum q>1): max |R*_4 - 12-term| = {d:.2e}')
# (d),(e) K=5 MILPs with structural restrictions
REPS=[m for m in range(32) if not m&1]
def milp5(q,b,zero_edges=()):
    nc=16; nb=16; N=nc+nb; R=[];lo=[]
    for t,m in enumerate(REPS):
        a=np.zeros(N); a[:nc]=F5[m]; a[nc+t]=-1; R.append(a); lo.append(0)
        a=np.zeros(N); a[:nc]=F5[31-m]; a[nc+t]=1; R.append(a); lo.append(1)
    lb=np.r_[-np.inf*np.ones(6),np.zeros(10),np.zeros(nb)]; ub=np.r_[np.inf*np.ones(nc),np.ones(nb)]
    for e in zero_edges: ub[6+e]=0
    r=milp(np.r_[1,q,b,np.zeros(nb)],constraints=LinearConstraint(np.array(R),lo,np.inf),integrality=np.r_[np.zeros(nc),np.ones(nb)],bounds=Bounds(lb,ub),options={'mip_rel_gap':1e-10}); return r.fun
iso=[e for e,(i,j) in enumerate(E5) if 0 in (i,j)]; rest=[e for e,(i,j) in enumerate(E5) if 0 not in (i,j)]
dd=0; nd=0; viol=0; ne=0
while time.time()-T0<230:
    q,b=sample(5,Z5,E5)
    lhs=milp5(q,b,iso); q4=q[1:]; b4=b[rest]; rhs=min(q[0],form4(q4,b4)); dd=max(dd,abs(lhs-rhs)); nd+=1
    S=set(rng.choice(5,rng.integers(1,3),replace=False)); within=[e for e,(i,j) in enumerate(E5) if (i in S)==(j in S)]
    v=milp5(q,b,within); ne+=1; viol+= v<q.min()-1e-7
print(f'(d) isolated-vertex class: {nd} K=5 instances, max |MILP - min(q_0,R*_4)| = {dd:.2e}')
print(f'(e) bipartite-support class: {ne} instances, MILP value < q_min in {viol}')
# (f)
import networkx as nx
from networkx.algorithms import isomorphism
bip=0; s3=0; s4=0; s5=0; reps=[]
for mask in range(1024):
    g=nx.Graph(); g.add_nodes_from(range(5)); g.add_edges_from(E5[k] for k in range(10) if mask>>k&1)
    if nx.is_bipartite(g): bip+=1; continue
    cov=len({v for e in g.edges for v in e})
    if cov==3: s3+=1
    elif cov==4: s4+=1
    else:
        s5+=1
        if not any(nx.is_isomorphic(g,h) for h in reps): reps.append(g)
print(f'(f) bipartite {bip}, non-bipartite covering 3/4/5 vertices: {s3}/{s4}/{s5}; unlabeled spanning non-bipartite: {len(reps)}')
# (g)
p=Fr(2,5); pi={0b10001:p}
for a,bb in itertools.combinations((1,2,3),2): pi[(1<<a)|(1<<bb)]=(1-p)/3
q=[sum(v for m,v in pi.items() if m>>i&1) for i in range(5)]
b=[sum(v for m,v in pi.items() if m>>i&1 and m>>j&1) for i,j in E5]
H=1-sum(q[:4])+sum(b[E5.index(e)] for e in itertools.combinations((1,2,3),2))+2*sum(b[E5.index((0,j))] for j in (1,2,3))
tri=[sum(b[E5.index(e)] for e in itertools.combinations(S,2)) for S in itertools.combinations(range(5),3)]
pent=[]
for r in itertools.permutations(range(1,5)):
    if r[0]>r[-1]: continue
    cyc=(0,)+r; pent.append(sum(b[E5.index(tuple(sorted((cyc[i],cyc[(i+1)%5]))))] for i in range(5)))
print(f'(g) p=2/5: q={[str(x) for x in q]}, sum={sum(q)}, sum-3qmin={sum(q)-3*min(q)}, H_T,0={H}, min odd cycle={min(tri+pent)}, MILP R*_5={milp5(np.array([float(x) for x in q]),np.array([float(x) for x in b])):.2e}')
print(f'[{time.time()-T0:.0f}s]')
