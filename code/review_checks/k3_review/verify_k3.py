"""Exact verification of D2_K3_27_orbit_Farkas_proof.json (format inferred from the K=4 verifier in [S17] App. B).
LP variables z=(pi[8], b[3], t). Rows (G z <= h), same order as the K=4 code:
  per edge k=(i,j): 3k: b<=q_i, 3k+1: b<=q_j, 3k+2: E[XiXj]<=b ; 9..11: q_i<=1/2 ; 12..15: t <= L_C - L_V
  for C in (dict 0, dict 1, dict 2, triangle).  Certificate: lambda>=0, mu with
  -e_t + G^T lambda + mu*e >= 0 on pi,b coords, =0 on t; bound=-lambda.h-mu >=0  =>  max t <= 0.
Then: (i) each Q is a vertex of D_F for some signed threshold rule F; (ii) every vertex of every rule
(exact cddlib) is an S3-image of a certified Q -> R*_3 = min{q1,q2,q3,b12+b13+b23} on the whole domain;
(iii) brute-force sanity check on random instances."""
import json, itertools, numpy as np
from fractions import Fraction as F
import cdd, cdd.gmp as G
from scipy.optimize import linprog
K=3; E=list(itertools.combinations(range(K),2)); X=[[(m>>i)&1 for i in range(K)] for m in range(8)]; NV=8+3+1
A=[[1,*x,*(x[i]*x[j] for i,j in E)] for x in X]+[[0]*4+[int(k==j) for k in range(3)] for j in range(3)]
C=[]
for i in range(3): c=[0]*7; c[1+i]=1; C.append(c)
C.append([0,0,0,0,1,1,1])
Gm=[];h=[]
for k,(i,j) in enumerate(E):
    for v in (i,j):
        z=[0]*NV
        for m in range(8): z[m]=-X[m][v]
        z[8+k]=1; Gm.append(z); h.append(F(0))
    z=[0]*NV
    for m in range(8): z[m]=X[m][i]*X[m][j]
    z[8+k]=-1; Gm.append(z); h.append(F(0))
for i in range(3):
    z=[0]*NV
    for m in range(8): z[m]=X[m][i]
    Gm.append(z); h.append(F(1,2))
REC=json.load(open('k3.json')); certified=[]
for n,rec in enumerate(REC):
    v=[F(s) for s in rec['Q']]; rows=[list(map(F,z)) for z in Gm]; hs=h.copy()
    for c in C:
        d=[F(u)-zz for u,zz in zip(c,v)]; z=[F(0)]*NV
        for m in range(8): z[m]=-sum(X[m][i]*d[1+i] for i in range(3))
        for k in range(3): z[8+k]=-d[4+k]
        z[11]=F(1); rows.append(z); hs.append(d[0])
    lam=[F(0)]*16
    for i,s in rec['dual_nonzero']: lam[i]=F(s)
    mu=F(rec['nu']); assert all(x>=0 for x in lam)
    red=[F(-1) if j==11 else F(0) for j in range(NV)]
    for i,l in enumerate(lam):
        if l:
            for j in range(NV): red[j]+=l*rows[i][j]
    for j in range(8): red[j]+=mu
    assert all(x>=0 for x in red[:11]) and red[11]==0, n
    bound=-sum(a*b for a,b in zip(lam,hs))-mu
    assert bound==F(rec['objective']) and bound>=0, n
    certified.append(tuple(v))
print(f'(0) {len(REC)} records: all Farkas certificates verify exactly (bound >= 0 => Delta <= 0); distinct Q: {len(set(certified))}')
# (i)/(ii) rules and vertices
def rule(w): return tuple(int(sum(w[i]*(2*x[i]-1) for i in range(3))>0) for x in X)
RULES=set()
for w0 in [(1,0,0),(1,1,1)]:
    for w in set(itertools.permutations(w0)):
        for sg in itertools.product((1,-1),repeat=3): RULES.add(rule(tuple(w[i]*sg[i] for i in range(3))))
print(f'    signed threshold rules for K=3: {len(RULES)}')
def perm_coef(c,p):
    beta=[0]*3; gam=[0]*3
    for i in range(3): beta[p[i]]=c[1+i]
    for e,(i,j) in enumerate(E): gam[E.index(tuple(sorted((p[i],p[j]))))]=c[4+e]
    return tuple([c[0]]+beta+gam)
cert_orbits={perm_coef(c,p) for c in certified for p in itertools.permutations(range(3))}
allv=set(); is_vertex_of_some_rule=set()
for f in RULES:
    rows=[[-f[m]]+A[m] for m in range(8)]+[[0]+A[8+j] for j in range(3)]
    g=G.copy_generators(G.polyhedron_from_matrix(G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY)))
    allv|={tuple(F(x) for x in r[1:]) for r in g.array if r[0]==1}
print(f'(i) every certified Q is a vertex of some rule polyhedron: {set(certified)<=allv}')
print(f'(ii) distinct vertices over all {len(RULES)} rules: {len(allv)}; all covered by S3-images of certified Q: {allv<=cert_orbits}; '
      f'S3-orbits of vertices: {len({min(perm_coef(v,p) for p in itertools.permutations(range(3))) for v in allv})}')
# (iii)
rng=np.random.default_rng(1); Z=np.array(X); P=np.array([Z[:,i]*Z[:,j] for i,j in E]); RL=[np.array(f,float) for f in RULES]; mx=0; nb=0
for _ in range(1500):
    k=rng.integers(2,8); S=rng.choice(8,k,replace=False); pi=np.zeros(8); pi[S]=rng.dirichlet(np.ones(k)); q=Z.T@pi
    if q.max()>=.5: continue
    J=P@pi; cap=np.array([min(q[i],q[j]) for i,j in E]); b=J+(cap-J)*rng.uniform(0,1,3)**2*(rng.uniform(size=3)<.5)
    tr=min(-linprog(-f,A_ub=P,b_ub=b,A_eq=np.vstack([np.ones(8),Z.T]),b_eq=np.r_[1,q],bounds=(0,None),method='highs').fun for f in RL)
    mx=max(mx,abs(tr-min(q.min(),b.sum()))); nb+=q.sum()>1
print(f'(iii) brute force (min over {len(RL)} rules of worst-case LP) vs min(q_min, b12+b13+b23): max diff {mx:.2e}; {nb} instances with sum q > 1')
# which certified Q are not threshold-rule vertices? test against all 16 self-dual rules
SD=[]
for bits in itertools.product((0,1),repeat=4):
    f=[0]*8
    for k,m in enumerate([m for m in range(8) if not m&1]): f[m]=bits[k]; f[7-m]=1-bits[k]
    SD.append(tuple(f))
vsd=set()
for f in SD:
    rows=[[-f[m]]+A[m] for m in range(8)]+[[0]+A[8+j] for j in range(3)]
    g=G.copy_generators(G.polyhedron_from_matrix(G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY)))
    vsd|={tuple(F(x) for x in r[1:]) for r in g.array if r[0]==1}
extra=[c for c in certified if c not in allv]
print(f'(iv) {len(extra)} certified Q are not threshold-rule vertices; of these, vertices of a (non-threshold) self-dual rule: {sum(c in vsd for c in extra)}; '
      f'valid certificates (Q>=0, complementary max>=1): {sum(all(sum(a*b for a,b in zip(c,A[m]))>=0 for m in range(8)) and all(max(sum(a*b for a,b in zip(c,A[m])),sum(a*b for a,b in zip(c,A[7-m])))>=1 for m in range(8)) for c in extra)}')
