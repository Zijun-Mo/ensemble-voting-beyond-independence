"""Independent completeness check of U248 on closure(R), bypassing the report's Lemma 2:
by threshold sufficiency (doc result 2(b)) every legal deterministic certificate lies in P_F for one of the 1882 signed
self-dual threshold rules F; up to S5 these are 60 rules whose P_F vertices (17,508, exact cdd) fall into 1660 S5-orbits.
For every orbit not in U248 we solve   max t  s.t.  t <= L_j(q(pi),b) - L_v(q(pi),b)  (j in U248),
pi in simplex, p_e(pi) <= b_e <= min(q_i,q_j), q_i <= 1/2, sum q - 3 q_k <= 1, b >= 0,
and convert the LP dual into an exact rational Farkas certificate of max t <= 0 (checked in Fractions).
Usage: python3 dominance_all.py <dir with report module> <start> <end>"""
import json, itertools, sys, time
import numpy as np
from fractions import Fraction as F
from scipy.optimize import linprog
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U
st,en=int(sys.argv[2]),int(sys.argv[3])
E=list(itertools.combinations(range(5),2)); Z=[[(m>>i)&1 for i in range(5)] for m in range(32)]
PAIR=[[Z[m][i]*Z[m][j] for (i,j) in E] for m in range(32)]
C=[c.coefficients() for c in U.make_u248()]
# variables: pi(32), b(10), t ; rows A z <= r
rows=[];rhs=[]
def row(pi,b,t,r): rows.append([F(x) for x in pi]+[F(x) for x in b]+[F(t)]); rhs.append(F(r))
for e,(i,j) in enumerate(E):
    z=[0]*10; z[e]=-1; row([PAIR[m][e] for m in range(32)],z,0,0)          # p_e - b_e <= 0
    for k in (i,j):
        z=[0]*10; z[e]=1; row([-Z[m][k] for m in range(32)],z,0,0)         # b_e - q_k <= 0
for k in range(5):
    row([Z[m][k] for m in range(32)],[0]*10,0,F(1,2))                       # q_k <= 1/2
    pass  # region R dropped: full feasible domain q_i <= 1/2
nbase=len(rows)
def Lrow(c):   # L_c as linear form in (pi,b) plus constant
    return [c[0]+sum(c[1+i]*Z[m][i] for i in range(5)) for m in range(32)], list(c[6:])
LC=[Lrow(c) for c in C]
orbs=[tuple(F(x) for x in o) for o in json.load(open('non_u248_orbits.json'))][st:en]
T=time.time(); worst=-1e9; bad=[]; nexact=0; maxden=1
for v in orbs:
    Lv=Lrow(v)
    R=[r[:] for r in rows]; RH=rhs[:]
    for (pj,bj) in LC:   # t - (L_j - L_v) <= 0, written with pi-coefficients (constants absorbed via sum pi = 1)
        R.append([-(pj[m]-Lv[0][m]) for m in range(32)]+[-(bj[e]-Lv[1][e]) for e in range(10)]+[F(1)]); RH.append(F(0))
    A=np.array([[float(x) for x in r] for r in R]); bvec=np.array([float(x) for x in RH])
    res=linprog(np.r_[np.zeros(42),-1.0],A_ub=A,b_ub=bvec,A_eq=np.r_[np.ones(32),np.zeros(11)][None,:],b_eq=[1.0],
                bounds=[(0,None)]*42+[(None,None)],method='highs')
    assert res.status==0
    tmax=-res.fun; worst=max(worst,tmax)
    mu=[F(float(-x)).limit_denominator(10**6) if x<-1e-12 else F(0) for x in res.ineqlin.marginals]
    nu=F(float(-res.eqlin.marginals[0])).limit_denominator(10**6)
    acc=[sum(mu[i]*R[i][k] for i in range(len(R)) if mu[i]) for k in range(43)]
    acc=[a+(nu if k<32 else 0) for k,a in enumerate(acc)]
    cost=sum(mu[i]*RH[i] for i in range(len(R)) if mu[i])+nu
    ok = acc[42]==1 and min(acc[:42])>=0 and cost<=0
    if ok: nexact+=1; maxden=max([maxden]+[x.denominator for x in mu if x]+[nu.denominator])
    else: bad.append((tmax,[str(x) for x in v]))
pos=[b for b in bad if b[0]>1e-9]
print('orbits with LP value t > 1e-9 (i.e. NOT dominated by U248 on the full domain):',len(pos), ' max t =',max([b[0] for b in bad],default=0))
print(f'orbits {st}..{st+len(orbs)}: max LP value of t = {worst:.3e}; exact Farkas certificates verified: {nexact}/{len(orbs)}; max denominator {maxden}; failures {len(bad)}  [{time.time()-T:.0f}s]')
for b in bad[:5]: print('  FAIL', b)
