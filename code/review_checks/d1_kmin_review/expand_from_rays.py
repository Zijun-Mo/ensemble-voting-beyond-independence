"""Given exact extreme rays (Normaliz, integer) of the tangent cone of Q(K_8) at a vertex v, compute every neighbour
v + t* r by an exact rational ratio test over all 3,636 rows, and classify each neighbour up to S_8 against the 23 orbits
found by the review's adjacency decomposition (run/adj_state_K8.json). Canonical form: refinement-based (from the report's
fast_canonical_rep8, an exact canonical form: it minimises over all permutations compatible with an isomorphism-invariant
vertex partition), cross-checked against the full 40320-permutation canonical form on a random sample.
Usage: python3 expand_from_rays.py <rays.json> <vertex as 'a|b|...' or 'F01'> <adj_state_K8.json> <report script>"""
import sys, json, itertools, math, time, ast, random
import numpy as np
from fractions import Fraction as F
rays=np.array(json.load(open(sys.argv[1])),dtype=np.int64)
n=8;E=list(itertools.combinations(range(n),2));idx={e:k for k,e in enumerate(E)}
rows=[]; rhs=[]
for L in (3,5,7):
    for S in itertools.combinations(range(n),L):
        for p in itertools.permutations(S[1:]):
            if p[0]>p[-1]: continue
            c=(S[0],)+p; r=np.zeros(28,dtype=np.int64)
            for i in range(L): r[idx[tuple(sorted((c[i],c[(i+1)%L])))]]=1
            rows.append(r); rhs.append(1)
for k in range(28):
    r=np.zeros(28,dtype=np.int64); r[k]=1; rows.append(r); rhs.append(0)
A=np.array(rows); b=np.array(rhs)
assert len(A)==3636
v=[F(1,3)]*28 if sys.argv[2]=='F01' else [F(x) for x in sys.argv[2].split('|')]
Lv=1
for x in v: Lv=Lv*x.denominator//math.gcd(Lv,x.denominator)
V=np.array([int(x*Lv) for x in v],dtype=np.int64)
slack=A@V-Lv*b; assert slack.min()>=0
tight=slack==0
# report's refinement canonical form (copied function; read before use)
src=open(sys.argv[4]).read(); s0=src.index('def fast_canonical_rep8'); s1=src.index('def exact_dd_fast128')
ns={'permutations':itertools.permutations,'EDGES8':E}; exec(src[s0:s1],ns); fcan=ns['fast_canonical_rep8']
perms=np.array(list(itertools.permutations(range(n))))
Pm=np.array([[idx[tuple(sorted((p[i],p[j])))] for (i,j) in E] for p in perms]);Pinv=np.argsort(Pm,axis=1)
def fullcanon(w):
    L=1
    for x in w: L=L*x.denominator//math.gcd(L,x.denominator)
    iv=np.array([int(x*L) for x in w]);W=iv[Pinv];c=W[np.lexsort(W.T[::-1])[0]]
    return tuple(F(int(x),L) for x in c)
St=json.load(open(sys.argv[3]))['orbits']
known={fcan(tuple(F(x) for x in k.split('|'))):k for k in St}
assert len(known)==len(St)
T=time.time(); AR=A@rays.T          # 3636 x nrays
nb=0; unb=0; new=set(); hits={}; sample=[]
for j in range(rays.shape[0]):
    ar=AR[:,j]
    assert np.all(ar[tight]>=0)
    neg=np.where(ar<0)[0]
    if len(neg)==0: unb+=1; continue
    s=slack[neg]; q=-ar[neg]; i=np.argmin(s/q); t=F(int(s[i]),int(q[i]))
    assert all(int(s[k])*int(q[i])>=int(s[i])*int(q[k]) for k in range(len(neg)))
    w=tuple(F(int(V[e]),Lv)+t*F(int(rays[j,e]),Lv) for e in range(28))
    c=fcan(w); nb+=1
    if c in known: hits[known[c]]=hits.get(known[c],0)+1
    else: new.add(c)
    if random.random()<0.002: sample.append((w,c))
for w,c in sample: assert (fullcanon(w)==fullcanon(c))
print(f'rays {rays.shape[0]}: unbounded {unb}, bounded edges {nb}; neighbour orbits hit {len(hits)}; NEW orbits {len(new)}  [{time.time()-T:.0f}s]; sample full-canon cross-checks {len(sample)} ok')
for c in list(new)[:5]: print('  NEW:',[str(x) for x in c])
