import json,pickle,itertools,numpy as np,math,sys
from fractions import Fraction as F
n=8;E=list(itertools.combinations(range(n),2));idx={e:k for k,e in enumerate(E)}
perms=np.array(list(itertools.permutations(range(n))))
P=np.array([[idx[tuple(sorted((p[i],p[j])))] for (i,j) in E] for p in perms]);Pinv=np.argsort(P,axis=1)
def canon(v):
    L=1
    for x in v: L=L*x.denominator//math.gcd(L,x.denominator)
    iv=np.array([int(x*L) for x in v]);W=iv[Pinv];w=W[np.lexsort(W.T[::-1])[0]]
    return '|'.join(str(F(int(x),L)) for x in w)
u7,uA,uB=pickle.load(open('c5k2.pkl','rb'))
S=json.load(open('run/adj_state_K8.json'))['orbits']
cA=canon([uA[e] for e in E]);cB=canon([uB[e] for e in E])
print('LiftA found:',cA in S,' LiftB found:',cB in S)
for k in S:
    v=[F(x) for x in k.split('|')]
    if all(x.denominator==1 for x in v): continue
    den=max(x.denominator for x in v)
    if den>3:
        tag='LiftA' if k==cA else 'LiftB' if k==cB else 'NEW'
        z=[E[i] for i,x in enumerate(v) if x==0]
        print(tag,'den',den,'zero edges',z,'\n   u =',{E[i]:str(x) for i,x in enumerate(v)})
