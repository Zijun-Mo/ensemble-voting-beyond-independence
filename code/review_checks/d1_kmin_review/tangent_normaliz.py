"""Extreme rays (Normaliz, exact) of the tangent cone of Q(K_8) at a given vertex. Usage: python3 tangent_normaliz.py <vertex 'a|b|..' or K44> <out.json>"""
import itertools, time, json, sys
import numpy as np
from fractions import Fraction as F
from PyNormaliz import Cone
n=8; E=list(itertools.combinations(range(n),2)); idx={e:k for k,e in enumerate(E)}
if sys.argv[1]=='K44': v=[F(1) if (a<4)==(b<4) else F(0) for a,b in E]
else: v=[F(x) for x in sys.argv[1].split('|')]
rows=[]
for L in (3,5,7):
    for S in itertools.combinations(range(n),L):
        for p in itertools.permutations(S[1:]):
            if p[0]>p[-1]: continue
            c=(S[0],)+p; es=[idx[tuple(sorted((c[i],c[(i+1)%L])))] for i in range(L)]
            if sum(v[e] for e in es)==1: rows.append([1 if k in es else 0 for k in range(28)])
for k in range(28):
    if v[k]==0: rows.append([1 if j==k else 0 for j in range(28)])
print('tight rows',len(rows),flush=True)
T=time.time(); R=Cone(inequalities=rows).ExtremeRays()
print('extreme rays:',len(R),f'[{time.time()-T:.1f}s]',flush=True)
json.dump(R,open(sys.argv[2],'w'))
