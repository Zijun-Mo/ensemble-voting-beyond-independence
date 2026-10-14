"""D2 fan graph: exact 8-variable LP for common symmetric adversary mass.

Graph E=(01,02,03,04,12,14,23); clique tree (014)-(012)-(023).
The algebraic projection is exact; scipy linprog is only a numerical solver.
Run: python D2_K5_fan_8var_exact_LP.py
"""
from itertools import combinations
import numpy as np
from scipy.optimize import linprog

E = [(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
TRI = [(0,1,2),(0,1,4),(0,2,3)]
STATES = np.array([[m >> i & 1 for i in range(5)] for m in range(32)], dtype=float)
PAIRS = np.array([STATES[:,i]*STATES[:,j] for i,j in E]).T

class Aff:
    """An affine function const + coef @ (t, u_01, ..., u_23)."""
    def __init__(self, const=0., coef=None):
        self.const = float(const)
        self.coef = np.zeros(8) if coef is None else np.array(coef, dtype=float)
    def __add__(self, other):
        if not isinstance(other, Aff): other=Aff(other)
        return Aff(self.const+other.const, self.coef+other.coef)
    __radd__=__add__
    def __neg__(self): return Aff(-self.const,-self.coef)
    def __sub__(self,other): return self+-other if isinstance(other,Aff) else self+(-other)
    def __rsub__(self,other): return -self+other


def fan_8var_lp(q, b):
    q=np.asarray(q,dtype=float); b=np.asarray(b,dtype=float)
    assert q.shape==(5,) and b.shape==(7,)
    v=[]
    for k in range(8):
        coef=np.zeros(8);coef[k]=1
        v.append(Aff(coef=coef))
    t=v[0];u=v[1:]
    r=[Aff(q[i])-t for i in range(5)]
    h=Aff(1)-t-t
    d=[Aff(b[k])-u[k] for k in range(7)]
    rows=[];rhs=[]
    def le(f,g=0):
        if not isinstance(g,Aff):g=Aff(g)
        diff=f-g
        rows.append(diff.coef);rhs.append(-diff.const)
    for i in range(5):
        le(t,q[i]);le(r[i],h)
    le(t+t,1)
    for k in range(7):
        le(u[k],t);le(u[k],b[k])
    # Symmetric 3-bit mass: for each clique, agreement edge masses lie
    # in the simplex of 3 one-agreement cuts and the all-agreement cut.
    for c in TRI:
        ids=[E.index(pair) for pair in combinations(c,2)]
        le(t,sum(u[k] for k in ids))
        for k in ids:
            le(sum(u[j] for j in ids if j!=k),u[k]+t)
    # Residual local pair Frechet inequalities on leaf edges and central 12.
    for k in [2,3,4,5,6]:
        i,j=E[k]
        le(r[i]+r[j]-h,d[k])
    # Projected interval for v01, v02, exact after eliminating 3rd moments.
    L01=[Aff(0),r[0]+r[1]-h,r[0]+r[1]+r[4]-h-d[3]-d[5]]
    U01=[r[0],r[1],d[0]]
    L02=[Aff(0),r[0]+r[2]-h,r[0]+r[2]+r[3]-h-d[2]-d[6]]
    U02=[r[0],r[2],d[1]]
    for aa in L01:
        for bb in U01:le(aa,bb)
    for aa in L02:
        for bb in U02:le(aa,bb)
    for aa in L01:
        for bb in L02:le(aa+bb,r[0]+d[4])
    for aa in U01:
        for bb in U02:le(r[0]+r[1]+r[2]-h-d[4],aa+bb)
    c=np.zeros(8);c[0]=-1
    return linprog(c,A_ub=np.array(rows),b_ub=np.array(rhs),bounds=[(0,None)]*8,method='highs')


def fan_full_lp(q,b):
    """Independent 48-variable primal (32 pi states +16 complement masses)."""
    q=np.asarray(q,float);b=np.asarray(b,float)
    Aeq=np.zeros((6,48));Aeq[0,:32]=1;Aeq[1:,:32]=STATES.T
    Aub=np.zeros((39,48));Aub[:7,:32]=PAIRS.T
    for j in range(16):
        for m in (j,31-j):Aub[7+m,m]=-1;Aub[7+m,32+j]=1
    c=np.zeros(48);c[32:]=-1
    return linprog(c,A_ub=Aub,b_ub=np.r_[b,np.zeros(32)],A_eq=Aeq,b_eq=np.r_[1,q],bounds=(0,None),method='highs')

if __name__=='__main__':
    rng=np.random.default_rng(175)
    count=0; old_condition_failed=0
    for k in range(800):
        q=np.array([.3,.32,.46,.46,.3])+rng.normal(0,.015,5)
        q=np.clip(q,.16,.49)
        if q[[0,1,4]].sum()>1:continue
        pi=np.zeros(32)
        for mask,prob in [(0,1-q[[0,1,4]].sum()),(1,q[0]),(2,q[1]),(16,q[4])]:
            for x2 in range(2):
                for x3 in range(2):
                    pi[mask|(x2<<2)|(x3<<3)] += prob*(q[2] if x2 else 1-q[2])*(q[3] if x3 else 1-q[3])
        indep=np.prod(np.where(STATES,q,1-q),axis=1)
        pi=(1-rng.uniform(.01,.5))*pi
        pi+= (1-pi.sum())*indep
        actual=PAIRS.T@pi
        b=actual+(np.minimum(q[[i for i,j in E]],q[[j for i,j in E]])-actual)*rng.uniform(0,1,7)
        a=fan_8var_lp(q,b); f=fan_full_lp(q,b)
        assert a.success and f.success,(a.message,f.message)
        assert abs(a.fun-f.fun)<1e-8,(q,b,a.fun,f.fun)
        cycles=[[(0,1),(0,2),(1,2)],[(0,1),(0,4),(1,4)],[(0,2),(0,3),(2,3)],[(0,3),(0,4),(1,2),(1,4),(2,3)]]
        candidate=min(*q,*[sum(b[E.index(e)] for e in cyc) for cyc in cycles])
        surplus=max(sum(q[list(tri)])-1 for tri in TRI)
        if candidate<surplus-1e-8:old_condition_failed+=1
        count+=1
    print('VERIFIED',count,'feasible cases: 8-variable = 48-variable LP; old sufficient condition fails in',old_condition_failed)
    print('The formulation is a theorem of clique gluing + triangle projection, not an inference from tests.')
