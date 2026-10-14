"""K=6: B_graph on the exactly enumerated fractional vertex orbits.
K=7,8: random-objective LP vertex sampling of Q(K_n); for every distinct fractional vertex found, compute B_graph(u)
by the exact judicious-partition MILP. Counts distinct vertices by a permutation-invariant signature."""
import sys, itertools, json, time, numpy as np
from fractions import Fraction as F
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
def setup(n):
    E=list(itertools.combinations(range(n),2)); idx={e:k for k,e in enumerate(E)}; cyc=set()
    for L in range(3,n+1,2):
        for S in itertools.combinations(range(n),L):
            for p in itertools.permutations(S[1:]):
                if p[0]>p[-1]: continue
                c=(S[0],)+p; cyc.add(frozenset(idx[tuple(sorted((c[i],c[(i+1)%L])))] for i in range(L)))
    A=np.zeros((len(cyc),len(E)))
    for r,C in enumerate(cyc): A[r,list(C)]=1
    return E,A
def bgraph(n,E,cost):
    m=len(E); reps=[S for S in range(1<<n) if not S&1]; nz=len(reps); N=m+nz; R=[];lo=[]
    for t,S in enumerate(reps):
        a=np.zeros(N); d=np.zeros(N)
        for k,(u,v) in enumerate(E):
            if (S>>u)&1 and (S>>v)&1: a[k]=1
            if not (S>>u)&1 and not (S>>v)&1: d[k]=1
        a[m+t]=-1; R.append(a); lo.append(0); d[m+t]=1; R.append(d); lo.append(1)
    r=milp(np.r_[cost,np.zeros(nz)],constraints=LinearConstraint(np.array(R),lo,np.inf),integrality=np.r_[np.zeros(m),np.ones(nz)],bounds=Bounds(0,1),options={'mip_rel_gap':1e-10,'time_limit':20})
    return r.fun
def signature(n,E,u):
    vals=tuple(sorted(np.round(u,6)))
    deg=tuple(sorted(tuple(sorted(np.round([u[k] for k,e in enumerate(E) if v in e],6))) for v in range(n)))
    return (vals,deg)
T0=time.time()
# K=6 exact orbits
E6,_=setup(6); orbs=json.load(open('Qfrac_orbits_K6.json'))
for o in orbs:
    u=np.array([float(F(x)) for x in o]); print('K=6 orbit: values',sorted(set(o),key=lambda s:F(s)),' zero edges',sum(F(x)==0 for x in o),' B_graph(u) =',round(bgraph(6,E6,u),6))
rng=np.random.default_rng(2)
for n,budget in ((7,120),(8,140)):
    E,A=setup(n); seen={}; nlp=0; T=time.time()
    while time.time()-T<budget:
        kind=rng.integers(3)
        c=rng.exponential(1,len(E)) if kind==0 else (rng.integers(1,6,len(E)).astype(float) if kind==1 else rng.uniform(0,1,len(E))**3+1e-3)
        r=linprog(c,A_ub=-A,b_ub=-np.ones(len(A)),bounds=(0,None),method='highs-ds'); nlp+=1
        u=r.x
        if np.all(np.abs(u*15-np.round(u*15))<1e-7) and np.all(np.abs(u-np.round(u))<1e-9): continue   # integral
        s=signature(n,E,u)
        if s not in seen: seen[s]=u
    res=[]
    for s,u in seen.items():
        if time.time()-T0>285: break
        b=bgraph(n,E,u); den=max((F(x).limit_denominator(60)).denominator for x in u); res.append((b,den,int(np.sum(u<1e-9))))
    bad=[r for r in res if r[0]<1-1e-6]
    print(f'K={n}: {nlp} random LPs, {len(seen)} distinct fractional vertex signatures (B computed for {len(res)}); '
          f'denominators seen {sorted(set(r[1] for r in res))}; zero-edge counts {sorted(set(r[2] for r in res))}; vertices with B_graph < 1: {len(bad)}')
print(f'[{time.time()-T0:.0f}s]')
