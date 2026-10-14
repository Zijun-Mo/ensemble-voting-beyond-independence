"""For every fractional vertex orbit found so far in the (partial) adjacency enumeration of Q(K_8):
matching-validity, smallest positive value, identification with the user's three vectors and with (1/3)b*(8,m),
and B_graph via the full-partition MILP (HiGHS; primal value and dual bound reported).
Also checks the user's three vectors exactly (Fractions): feasibility and rank of tight rows (sympy)."""
import json,pickle,itertools,math,sys,numpy as np,sympy as sp
from fractions import Fraction as F
from scipy.optimize import milp, LinearConstraint, Bounds
n=8;E=list(itertools.combinations(range(n),2));m=len(E);idx={e:k for k,e in enumerate(E)}
perms=np.array(list(itertools.permutations(range(n))))
P=np.array([[idx[tuple(sorted((p[i],p[j])))] for (i,j) in E] for p in perms]);Pinv=np.argsort(P,axis=1)
def canon(v):
    L=1
    for x in v: L=L*x.denominator//math.gcd(L,x.denominator)
    iv=np.array([int(x*L) for x in v]);W=iv[Pinv];w=W[np.lexsort(W.T[::-1])[0]]
    return '|'.join(str(F(int(x),L)) for x in w)
cyc=set()
for L in (3,5,7):
    for S in itertools.combinations(range(n),L):
        for p in itertools.permutations(S[1:]):
            if p[0]>p[-1]: continue
            c=(S[0],)+p; cyc.add(frozenset(idx[tuple(sorted((c[i],c[(i+1)%L])))] for i in range(L)))
cyc=[sorted(C) for C in cyc]
def exact_vertex(v):
    sl=[sum(v[k] for k in C)-1 for C in cyc]
    if min(sl)<0 or min(v)<0: return False,None
    rows=[[1 if k in C else 0 for k in range(m)] for C,s in zip(cyc,sl) if s==0]+[[1 if k==j else 0 for k in range(m)] for j in range(m) if v[j]==0]
    return True, sp.Matrix(rows).rank()
reps=[S for S in range(1<<n) if not S&1]
Rm=[];lo=[]
for t,S in enumerate(reps):
    a=np.zeros(m+len(reps)); d=np.zeros(m+len(reps))
    for k,(x,y) in enumerate(E):
        if (S>>x)&1 and (S>>y)&1: a[k]=1
        if not (S>>x)&1 and not (S>>y)&1: d[k]=1
    a[m+t]=-1; Rm.append(a); lo.append(0); d[m+t]=1; Rm.append(d); lo.append(1)
Rm=np.array(Rm)
def bgraph(v):
    r=milp(np.r_[np.array([float(x) for x in v]),np.zeros(len(reps))],constraints=LinearConstraint(Rm,lo,np.inf),
           integrality=np.r_[np.zeros(m),np.ones(len(reps))],bounds=Bounds(0,1),options={'mip_rel_gap':1e-12})
    return r.fun,r.mip_dual_bound
def bstar(mz):
    bot={2*i+1 for i in range(mz)}; z={(2*i,2*i+1) for i in range(mz)}
    return [F(0) if e in z else (F(1,3) if (e[0] in bot)==(e[1] in bot) else F(2,3)) for e in E]
canon_b={canon(bstar(mz)):mz for mz in range(4)}
U=pickle.load(open('user_k8/user_vecs.pkl','rb'))
ucanon={}
for k,v in U.items():
    ok,rk=exact_vertex(list(v)); ucanon[canon(list(v))]=k
    print(f'user vector {k}: exact feasible={ok}, exact tight rank={rk}/28')
S=json.load(open('run/adj_state_K8.json'))['orbits']
print('orbits found so far:',len(S),' unexpanded:',sum(not x['expanded'] for x in S.values()))
for key,info in S.items():
    v=[F(x) for x in key.split('|')]
    if all(x.denominator==1 for x in v): continue
    z=[E[i] for i,x in enumerate(v) if x==0]; vs=[a for e in z for a in e]; match=len(vs)==len(set(vs))
    den=max(x.denominator for x in v); mn=min(x for x in v if x>0)
    tag=ucanon.get(key,'')+(' canonical b*(8,%d)/3'%canon_b[key] if key in canon_b else '')
    b,db=bgraph(v)
    print(f'den {den}  zero edges {len(z)} {"matching" if match else "NON-matching"}  min positive {mn}  expanded={info["expanded"]}  {tag:28s} B={b:.9f} dual={db:.9f}')
