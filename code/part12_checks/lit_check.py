import numpy as np
rng=np.random.default_rng(1)
def thm3(p):
    K=len(p); m=(K+1)//2; s=np.sort(p)[::-1]
    return max(0,max((s[:j].sum()-(m-1))/(j-m+1) for j in range(m,K+1)))
def matan_lower(p):  # p ascending, k=(L+1)/2, xi(m)=(1/m)sum_{i=k-m+1}^{L} p_i - (L-k)/m
    L=len(p); k=(L+1)//2; a=np.sort(p)
    return max([0]+[a[k-m:].sum()/m-(L-k)/m for m in range(1,k+1)])
def ruger_upper(alpha,k):  # max P(at least k of n) = min{1, min_{0<=j<=k-1} sum_{i=1}^{n-j} alpha_(i)/(k-j)}, ascending
    n=len(alpha); a=np.sort(alpha)
    return min([1]+[a[:n-j].sum()/(k-j) for j in range(k)])
e1=e2=0
for K in [1,3,5,7,9,11]:
    for _ in range(2000):
        p=rng.uniform(0,1,K)
        e1=max(e1,abs(thm3(p)-matan_lower(p)))
        e2=max(e2,abs(thm3(p)-(1-ruger_upper(1-p,(K+1)//2))))
print("Theorem 3 vs Matan lower bound, max diff:",e1)
print("Theorem 3 vs 1 - Rueger-type bound on errors, max diff:",e2)

# Kaniovski (2010) Table 2, n=3
def thm4(p,c):
    r=p*p+c*p*(1-p); return max(r,6*p-3*r-2), min(3*r,2*p-r)
def sec27(p,c):
    r={}
    for i,j in [(0,1),(0,2),(1,2)]: r[i,j]=p[i]*p[j]+c*np.sqrt(p[i]*(1-p[i])*p[j]*(1-p[j]))
    R=sum(r.values()); L=max(0,r[0,1]+r[0,2]-p[0],r[0,1]+r[1,2]-p[1],r[0,2]+r[1,2]-p[2])
    U=min(r[0,1],r[0,2],r[1,2],1-sum(p)+R); return R-2*U, R-2*L
for (p,c,tab) in [(.5,0,(.25,.75)),(.5,.2,(.3,.7)),(.5,-.2,(.4,.6)),(.75,0,(.813,.938)),(.75,.2,(.7,.9))]:
    lo,hi=thm4(p,c); print(f"p={p},c={c}: Thm4=[{lo:.4f},{hi:.4f}]  Kaniovski Table2={tab}")
lo,hi=sec27([.75,.6,.6],.2); print(f"p=(.75,.6,.6),c=.2: Sec2.7=[{lo:.4f},{hi:.4f}]  Kaniovski Table2=(0.577,0.792)")

# Zaigraev & Kaniovski (2012): min over all dependence = (n(2p-1)+1)/(n+1), attained at x2=(2p-1)/(2p); max=min(1,2np/(n+1))
for p in [0.6,0.75,0.9]:
    n=3; x2=(2*p-1)/(2*p); r=p*p+x2*p*(1-p)
    print(f"p={p}: ZK min={(n*(2*p-1)+1)/(n+1):.4f}, Thm4 lower at ZK's x2={thm4(p,x2)[0]:.4f}, r=(3p-1)/2? {abs(r-(3*p-1)/2)<1e-12}")
