"""Independent checks for the K5 odd-cycle / positive-lift report.
(1) exact vertex enumeration of the K5 odd-cycle covering polyhedron (cddlib, GMP rationals);
(2) Lemma 2.1: exact check of the 2+3 split bound on random rational y (and extreme y);
(3) Theorem 2.3: B_graph (pure-quadratic certificate MILP) == c* on random b;
(4) Theorem 3.1: R*_det (full MILP) == min{q_i, c*} on random R0 instances;
(5) Prop 5.1: admissible intercept set, exact rational check at all breakpoints and midpoints."""
import itertools, time, numpy as np
from fractions import Fraction as Fr
import cdd, cdd.gmp as G
from scipy.optimize import milp, LinearConstraint, Bounds
from exact_core import Z, FEAT, EDGES, REPS, rstar
T0=time.time(); rng=np.random.default_rng(3)
tri=[tuple(EDGES.index(e) for e in itertools.combinations(S,2)) for S in itertools.combinations(range(5),3)]
pent=[]
for r in itertools.permutations(range(1,5)):
    if r[0]>r[-1]: continue
    c=(0,)+r; pent.append(tuple(sorted(EDGES.index(tuple(sorted((c[i],c[(i+1)%5])))) for i in range(5))))
CYC=tri+pent; assert len(set(CYC))==22
# (1)
rows=[[-1]+[1 if e in C else 0 for e in range(10)] for C in CYC]+[[0]+[1 if k==e else 0 for k in range(10)] for e in range(10)]
g=G.copy_generators(G.polyhedron_from_matrix(G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY)))
V=[tuple(Fr(x) for x in r[1:]) for r in g.array if r[0]==1]
frac=[v for v in V if any(x.denominator>1 for x in v)]
print(f'(1) odd-cycle cover polyhedron of K5: {len(V)} vertices, fractional: {[tuple(str(x) for x in v) for v in frac]}')
# (2)
SUBS=[m for m in range(32) if bin(m).count('1')==2]
def inner(y,m): return sum(y[k] for k,(i,j) in enumerate(EDGES) if (m>>i)&1 and (m>>j)&1)
bad=0
for t in range(20000):
    mode=t%3
    y=[Fr(int(v)) for v in rng.integers(0,20,10)] if mode==0 else [Fr(int(v)) for v in rng.integers(0,3,10)*rng.integers(1,1000,10)]
    if mode==2: y=[Fr(0)]*10; 
    if mode==2:
        for k in rng.choice(10,rng.integers(1,4),replace=False): y[k]=Fr(int(rng.integers(1,100)))
    Y=sum(y)
    if Y==0: continue
    if not any(inner(y,m)<=Y/3 and inner(y,31^m)<=Y/3 for m in SUBS): bad+=1
print(f'(2) 2+3 split with both sides <= Y/3 failed in {bad} of 20000 exact random y')
# (3) B_graph via pure-quadratic certificate MILP vs c*
def bgraph(b):
    nc=16; nb=16; n=nc+nb; R=[];lo=[]
    for t,m in enumerate(REPS):
        a=np.zeros(n); a[:nc]=FEAT[m]; a[nc+t]=-1; R.append(a); lo.append(0)
        a=np.zeros(n); a[:nc]=FEAT[31-m]; a[nc+t]=1; R.append(a); lo.append(1)
    lb=np.zeros(n); ub=np.r_[np.zeros(6),np.inf*np.ones(10),np.ones(nb)]
    r=milp(np.r_[0,np.zeros(5),b,np.zeros(nb)],constraints=LinearConstraint(np.array(R),lo,np.inf),integrality=np.r_[np.zeros(nc),np.ones(nb)],bounds=Bounds(lb,ub),options={'mip_rel_gap':1e-11}); return r.fun
cst=lambda b: min(sum(b[e] for e in C) for C in CYC)
d3=0; n3=0
while time.time()-T0<120:
    b=rng.uniform(0,1,10)**rng.choice([1,3,6]); b*=(rng.uniform(size=10)<rng.choice([0.6,1.0]))
    d3=max(d3,abs(bgraph(b)-cst(b))); n3+=1
print(f'(3) B_graph vs c*: {n3} random b, max |diff| = {d3:.2e}')
# (4) full R*_det vs 27-term formula in R0
d4=0; n4=0
while time.time()-T0<220:
    k=rng.integers(2,12); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k)); q=Z.T@pi
    if q.min()<0.01 or q.max()>=0.5 or q.sum()>1: continue
    J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
    b=J+(cap-J)*rng.uniform(0,1,10)**2*(rng.uniform(size=10)<0.5)
    d4=max(d4,abs(rstar(q,b)[0]-min(q.min(),cst(b)))); n4+=1
print(f'(4) R0 instances: {n4}, max |R*_MILP - min(q_min,c*)| = {d4:.2e}')
# (5) admissible intercepts
def H(s,gm,m): return sum(s[i] for i in range(5) if (m>>i)&1)+inner(gm,m)
def legal(a,s,gm):
    if a<0: return False
    Q=[H(s,gm,m)-a*(bin(m).count('1')-1) for m in range(32)]
    return all(v>=0 for v in Q) and all(max(Q[m],Q[31^m])>=1 for m in range(16))
def formula(s,gm):
    sz=lambda m: bin(m).count('1')
    A=min([H(s,gm,m)/(sz(m)-1) for m in range(32) if sz(m)>=2]+[(H(s,gm,31^(1<<i))-1)/3 for i in range(5) if s[i]<1]
          +[max(H(s,gm,m)-1,(H(s,gm,31^m)-1)/2) for m in range(32) if sz(m)==2])
    T=(H(s,gm,31)-1)/4; return A,T
mism=0; tested=0
for _ in range(400):
    s=[Fr(int(v),4) for v in rng.integers(0,7,5)]; gm=[Fr(int(v),4) for v in rng.integers(0,7,10)]
    if not legal(Fr(0),s,gm): continue
    A,T=formula(s,gm); tested+=1
    pts=sorted({Fr(0),A,T,Fr(1),A+Fr(1,7),T+Fr(1,7),Fr(1)-Fr(1,7),Fr(1)+Fr(1,7)}|{x/2 for x in (A,T,Fr(1))})
    pts=sorted(set(pts)|{(u+v)/2 for u,v in zip(pts,pts[1:])})
    for a in pts:
        pred=(0<=a<=min(A,T)) or (A>=1 and 1<=a<=A)
        if pred!=legal(a,s,gm): mism+=1
print(f'(5) intercept formula: {tested} certificates in D0, mismatches at breakpoints/midpoints: {mism}   [{time.time()-T0:.0f}s]')
