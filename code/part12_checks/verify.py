import itertools, numpy as np
from math import comb
from scipy.optimize import linprog
rng = np.random.default_rng(0)

def states(K): return np.array(list(itertools.product([0,1], repeat=K)))

def lp_bound(p, r=None, sense='min'):
    K=len(p); Z=states(K); m=(K+1)//2
    c=(Z.sum(1)>=m).astype(float)
    A=[np.ones(len(Z))]; b=[1.0]
    for i in range(K): A.append(Z[:,i]); b.append(p[i])
    if r is not None:
        for i in range(K):
            for j in range(i+1,K): A.append(Z[:,i]*Z[:,j]); b.append(r[i][j])
    res=linprog(c if sense=='min' else -c, A_eq=np.array(A), b_eq=b, bounds=(0,None), method='highs')
    if not res.success: return None
    return res.fun if sense=='min' else -res.fun

# --- Binomial monotonicity formula
def Aind(K,p): return sum(comb(K,j)*p**j*(1-p)**(K-j) for j in range((K+1)//2,K+1))
err=0
for K in [1,3,5,7,9,11]:
    r=(K-1)//2
    for p in rng.uniform(0,1,20):
        lhs=Aind(K+2,p)-Aind(K,p); rhs=comb(K,r)*p**(r+1)*(1-p)**(r+1)*(2*p-1)
        err=max(err,abs(lhs-rhs))
print("binomial increment formula max err", err)

# --- Theorem 2 & 3
def thm3(p):
    K=len(p); m=(K+1)//2; s=np.sort(p)[::-1]
    return max(0,max((s[:j].sum()-(m-1))/(j-m+1) for j in range(m,K+1)))
err2=err3=0
for K in [1,3,5,7,9]:
    m=(K+1)//2
    for p0 in rng.uniform(0,1,15):
        f=max(0,(2*K*p0-K+1)/(K+1)); err2=max(err2,abs(lp_bound([p0]*K)-f))
    for _ in range(40):
        p=rng.uniform(0,1,K)
        if rng.random()<0.3: p=rng.uniform(0.4,1,K)
        err3=max(err3,abs(lp_bound(p)-thm3(p)))
print("Theorem 2 max err", err2); print("Theorem 3 max err", err3)

# --- Theorem 4 / feasibility (8) for K=3 symmetric
err4=0; feas_mismatch=0; n=0
for _ in range(3000):
    p=rng.uniform(0.01,0.99); rho=rng.uniform(-1,1)
    r=p*p+rho*p*(1-p); R=[[r]*3]*3
    lo=lp_bound([p]*3,R,'min')
    pred_feas = max(0,p-1/3,2*p-1)-1e-12<=r<=p+1e-12
    if (lo is not None)!=pred_feas: feas_mismatch+=1; continue
    if lo is None: continue
    hi=lp_bound([p]*3,R,'max'); n+=1
    err4=max(err4,abs(lo-max(r,6*p-3*r-2)),abs(hi-min(3*r,2*p-r)))
print(f"Theorem 4 max err {err4} over {n} feasible cases; feasibility-condition mismatches {feas_mismatch}")

# --- Section 2.7 general K=3
err7=0;n7=0;mis7=0
for _ in range(3000):
    p=rng.uniform(0.05,0.95,3); r=np.zeros((3,3))
    for i,j in [(0,1),(0,2),(1,2)]:
        rho=rng.uniform(-1,1); r[i][j]=p[i]*p[j]+rho*np.sqrt(p[i]*(1-p[i])*p[j]*(1-p[j]))
    r12,r13,r23=r[0][1],r[0][2],r[1][2]; Rs=r12+r13+r23
    L=max(0,r12+r13-p[0],r12+r23-p[1],r13+r23-p[2]); U=min(r12,r13,r23,1-p.sum()+Rs)
    lo=lp_bound(p,r,'min')
    if (lo is not None)!=(L<=U+1e-12): mis7+=1; continue
    if lo is None: continue
    hi=lp_bound(p,r,'max'); n7+=1
    err7=max(err7,abs(lo-(Rs-2*U)),abs(hi-(Rs-2*L)))
print(f"Section 2.7 max err {err7} over {n7} feasible cases; feasibility mismatches {mis7}")

# --- Section 2.6 example
p=0.6; r=0.36; R=[[r]*3]*3
print("p=0.6,rho=0 LP min/max:", lp_bound([p]*3,R,'min'), lp_bound([p]*3,R,'max'), " independent:", Aind(3,0.6))
