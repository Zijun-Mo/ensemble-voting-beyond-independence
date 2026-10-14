"""Extreme rays of the tangent cone of Q(K_8) at F01 = (1/3)*1, i.e. T = {d : d(triangle) >= 0 for all 56 triangles}, via Normaliz."""
import itertools, time, json, sys
from PyNormaliz import Cone
n=8; E=list(itertools.combinations(range(n),2)); idx={e:k for k,e in enumerate(E)}
rows=[[1 if k in (idx[(a,b)],idx[(a,c)],idx[(b,c)]) else 0 for k in range(28)] for a,b,c in itertools.combinations(range(n),3)]
T=time.time()
C=Cone(inequalities=rows)
R=C.ExtremeRays()
print('extreme rays:',len(R),f'[{time.time()-T:.1f}s]',flush=True)
json.dump(R,open(sys.argv[1],'w'))
