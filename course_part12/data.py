import itertools, numpy as np
from math import comb
from scipy.optimize import linprog
rng=np.random.default_rng(7)
def Aind(K,p): return sum(comb(K,j)*p**j*(1-p)**(K-j) for j in range((K+1)//2,K+1))
# Part I curve data p=0.7
p=0.7
print("% K, Aind, Aminus, Aplus")
for K in range(1,22,2):
    print(f"({K},{Aind(K,p):.4f}) ({K},{max(0,(2*K*p-K+1)/(K+1)):.4f}) ({K},{min(1,2*K*p/(K+1)):.4f})")

def states(K): return np.array(list(itertools.product([0,1],repeat=K)))
def full_lp(p,r,sense):
    K=len(p); Z=states(K); m=(K+1)//2; c=(Z.sum(1)>=m).astype(float)
    A=[np.ones(len(Z))]; b=[1.0]
    for i in range(K): A.append(Z[:,i]); b.append(p[i])
    if r is not None:
        for i in range(K):
            for j in range(i+1,K): A.append(Z[:,i]*Z[:,j]); b.append(r[i][j])
    res=linprog(c if sense=='min' else -c,A_eq=np.array(A),b_eq=b,bounds=(0,None),method='highs')
    return None if not res.success else (res.fun if sense=='min' else -res.fun)
def reduced_lp(K,p,r,sense):
    s=np.arange(K+1); m=(K+1)//2; c=(s>=m).astype(float)
    A=np.vstack([np.ones(K+1), s, s*(s-1)/2]); b=[1,K*p,comb(K,2)*r]
    res=linprog(c if sense=='min' else -c,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
    if not res.success: return None, None
    return (res.fun if sense=='min' else -res.fun), res.x

# (a) upper bound for unequal p (Matan upper) vs LP
def upper(p):
    K=len(p); m=(K+1)//2; a=np.sort(p)
    return min([1]+[a[:j].sum()/(j-m+1) for j in range(m,K+1)])
e=0
for K in [1,3,5,7]:
    for _ in range(60):
        pp=rng.uniform(0,1,K); e=max(e,abs(full_lp(pp,None,'max')-upper(pp)))
print("upper bound unequal p vs LP:",e)

# (b) systematic sampling construction check
def syst(x,U):
    c=np.concatenate([[0],np.cumsum(x)])
    return np.array([int(np.floor(c[i+1]-U)-np.floor(c[i]-U))  for i in range(len(x))])  # #points of U+Z in [c_i,c_{i+1})
e=0;bad=0
for _ in range(200):
    K=rng.integers(1,9); x=rng.uniform(0,1,K); Us=(np.arange(20000)+0.5)/20000
    Zs=np.array([syst(x,U) for U in Us]); e=max(e,np.abs(Zs.mean(0)-x).max())
    s=x.sum(); tot=Zs.sum(1); bad+=np.any((tot<np.floor(s))|(tot>np.ceil(s)))|np.any(Zs>1)
print("systematic sampling marginal err",e,"violations",bad)

# (c) reduced LP == full LP for symmetric K=5
e=0
for _ in range(40):
    pp=rng.uniform(0.05,0.95); rho=rng.uniform(-0.3,1); r=pp*pp+rho*pp*(1-pp)
    for sense in ['min','max']:
        f=full_lp([pp]*5,[[r]*5]*5,sense); g,_=reduced_lp(5,pp,r,sense)
        if (f is None)!=(g is None): print("feas mismatch"); continue
        if f is not None: e=max(e,abs(f-g))
print("reduced vs full LP, K=5:",e)

# (d) table for general K
for (pp,rho) in [(0.6,0.0),(0.6,0.2),(0.7,0.1)]:
    r=pp*pp+rho*pp*(1-pp); row=[]
    for K in [3,5,7,9,11,21]:
        lo,xl=reduced_lp(K,pp,r,'min'); hi,xh=reduced_lp(K,pp,r,'max')
        row.append(f"K={K}: [{lo:.4f}, {hi:.4f}] supp_min={np.nonzero(xl>1e-9)[0].tolist()} supp_max={np.nonzero(xh>1e-9)[0].tolist()}" if lo is not None else f"K={K}: infeasible")
    print(f"p={pp}, rho={rho}, ind(K)={[round(Aind(K,pp),4) for K in [3,5,7,9,11,21]] if rho==0 else ''}")
    for x in row: print("   ",x)
