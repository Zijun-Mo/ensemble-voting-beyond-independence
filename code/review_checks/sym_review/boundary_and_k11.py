"""Usage: python3 -I boundary_and_k11.py boundary|k11  (boundary: K=5..8 boundary points vs MILP; k11: K=11 vertices with the [S31] symmetry-broken MILP)"""
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE,'..','..','..','sources','S31_bundle'))
import time; mode=sys.argv[1]; sys.argv=[sys.argv[0],'none']
import check_s31 as C
from fractions import Fraction as Fr
import numpy as np
T0=time.time()
if mode=='boundary':
    for K in (5,6,7,8):
        E,Z,P=C.setup(K); Tr=C.reduced_family(K); w=0; n=0
        for i in range(1,40):
            q=Fr(i,80)
            for b in (C.bmin(K,q), q):              # lower and upper boundary
                v,_=C.det_milp(K,float(q),float(b),E,Z,P); w=max(w,abs(v-float(C.U(Tr,q,b)))); n+=1
        for j in range(0,21):                       # right boundary q = 1/2
            q=Fr(1,2); lo=C.bmin(K,q); b=lo+(q-lo)*Fr(j,20)
            v,_=C.det_milp(K,float(q),float(b),E,Z,P); w=max(w,abs(v-float(C.U(Tr,q,b)))); n+=1
        print(f'K={K}: {n} boundary points, max |R_det - U_det| = {w:.1e}  [{time.time()-T0:.0f}s]')
else:
    import verify_symmetric_closed_form as V   # the report's own symmetry-broken MILP
    K=11; Tf=C.full_family(K); Tr=C.reduced_family(K)
    vs=[v for v in C.cell_vertices(K,Tf) if v[2]>0 and v[0]<Fr(1,2)]
    print(f'K=11: {len(vs)} nontrivial cell vertices (independently enumerated)')
    for q,b,u,tight,onb in vs:
        val,gap=V.exact_milp_sym_fast(K,q,b,time_limit=60)
        print(f'   ({q},{b}) U={u} ({float(u):.9f})  MILP={val:.9f}  gap={gap}  tight={tight}  [{time.time()-T0:.0f}s]', flush=True)
        if time.time()-T0>240: print('   time limit'); break
