"""Exact orbit-wise vertex enumeration of Q(K_n)={u>=0: u(C)>=1 for all odd cycles C} by adjacency decomposition.
The vertex graph of a pointed polyhedron is connected, so starting from one vertex and repeatedly taking all
neighbours (extreme rays of the tangent cone, computed exactly by cddlib/GMP; exact integer ratio test) of one
representative per S_n-orbit reaches every vertex orbit.  Resumable: state in adj_state_K{n}.json.  Each run stops
after a time budget.  Usage: python3 adj_orbits.py n [budget_seconds]"""
import sys, itertools, time, json, math, os
from fractions import Fraction as F
import numpy as np
import cdd, cdd.gmp as G
n=int(sys.argv[1]); budget=float(sys.argv[2]) if len(sys.argv)>2 else 260
SKIP_INT=len(sys.argv)>3 and sys.argv[3]=='skipint'
T0=time.time()
E=list(itertools.combinations(range(n),2)); m=len(E); idx={e:k for k,e in enumerate(E)}
cyc=set()
for L in range(3,n+1,2):
    for S in itertools.combinations(range(n),L):
        for p in itertools.permutations(S[1:]):
            if p[0]>p[-1]: continue
            c=(S[0],)+p; cyc.add(frozenset(idx[tuple(sorted((c[i],c[(i+1)%L])))] for i in range(L)))
A=np.zeros((len(cyc),m),dtype=np.int64)
for r,C in enumerate(sorted(cyc,key=sorted)): A[r,list(C)]=1
# full system  A u >= 1 , I u >= 0
AA=np.vstack([A,np.eye(m,dtype=np.int64)]); rhs=np.r_[np.ones(len(A),dtype=np.int64),np.zeros(m,dtype=np.int64)]
perms=np.array(list(itertools.permutations(range(n))))
P=np.array([[idx[tuple(sorted((p[i],p[j])))] for (i,j) in E] for p in perms])   # edge k -> P[s,k]
Pinv=np.argsort(P,axis=1)
def to_int(v):
    L=1
    for x in v: L=L*x.denominator//math.gcd(L,x.denominator)
    return L,np.array([int(x*L) for x in v],dtype=np.int64)
def canon(v):
    L,iv=to_int(v); W=iv[Pinv]            # W[s,k] = v[Pinv[s,k]] : permuted vector
    order=np.lexsort(W.T[::-1]); w=W[order[0]]
    return tuple(str(F(int(x),L)) for x in w)
def neighbours(v):
    L,U=to_int(v); slack=AA@U-L*rhs
    assert slack.min()>=0
    I=np.where(slack==0)[0]
    rows=[[0]+[int(x) for x in AA[i]] for i in I]
    H=G.matrix_from_array(rows,rep_type=cdd.RepType.INEQUALITY)
    g=G.copy_generators(G.polyhedron_from_matrix(H))
    assert not g.lin_set, 'tangent cone not pointed'
    out=[]
    for r in g.array:
        if r[0]!=0: continue
        d=[F(x) for x in r[1:]]
        Ld=1
        for x in d: Ld=Ld*x.denominator//math.gcd(Ld,x.denominator)
        D=np.array([int(x*Ld) for x in d],dtype=np.int64)
        ad=AA@D; neg=np.where(ad<0)[0]
        if len(neg)==0: continue           # unbounded edge
        # t = slack_i/(L*(-ad_i)) (in units: v + t*D/Ld ... ) use exact Fractions on candidates
        s=slack[neg]; q=-ad[neg]
        j=np.argmin(s/q); sj,qj=int(s[j]),int(q[j])
        # exact check: s_i*qj >= sj*q_i for all i ; then minimisers
        lhs=s.astype(object)*qj; rhs_=q.astype(object)*sj
        assert all(lhs>=rhs_)
        t=F(sj,qj)   # in scaled units: U + t*D  (U = L*v)
        w=tuple(F(int(U[k]),L)+t*F(int(D[k]),L) for k in range(m))
        out.append(w)
    return out, len(I), len(g.array)
st=f'adj_state_K{n}.json'
if os.path.exists(st):
    S=json.load(open(st))
else:
    v0=tuple(F(0) if (a==0)!=(b==0) else F(1) for (a,b) in E)   # 1 - chi(delta({0}))
    S={'orbits':{ '|'.join(canon(v0)):{'expanded':False}} }
orb=S['orbits']
while True:
    todo=sorted([k for k,x in orb.items() if not x['expanded']],key=lambda k:(all(F(x).denominator==1 for x in k.split('|')),k.split('|').count('0')))
    if not todo: break
    if SKIP_INT:
        todo=[k for k in todo if not all(F(x).denominator==1 for x in k.split('|')) and len(set(k.split('|')))>1]
        if not todo: break
    if time.time()-T0>budget: break
    k=todo[0]; v=tuple(F(x) for x in k.split('|'))
    t=time.time(); nb,nt,nr=neighbours(v)
    new=0
    for w in nb:
        c='|'.join(canon(w))
        if c not in orb: orb[c]={'expanded':False}; new+=1
    orb[k].update(expanded=True,tight=nt,rays=nr)
    print(f'expanded orbit (tight {nt}, cone rays {nr}) -> {len(nb)} nbrs, {new} new orbits; total {len(orb)}  [{time.time()-t:.1f}s]',flush=True)
    json.dump(S,open(st,'w'))
json.dump(S,open(st,'w'))
todo=[k for k,x in orb.items() if not x['expanded']]
fr=[k for k in orb if any(F(x).denominator>1 for x in k.split('|'))]
print(f'K={n}: {len(orb)} vertex orbits found, {len(fr)} fractional; unexpanded {len(todo)}; COMPLETE={not todo}  [{time.time()-T0:.0f}s]')
for k in fr:
    vals=sorted(set(k.split('|')),key=F); z=sum(F(x)==0 for x in k.split('|'))
    print('  frac orbit: values',vals,' zero edges',z)
