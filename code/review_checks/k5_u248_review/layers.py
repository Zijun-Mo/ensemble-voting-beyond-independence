"""Report §4 table: number of irredundant U248 terms on the saturation layers D_S (b_e = p_e(pi) for e in S), within R.
Term j is necessary on a layer iff  max_{layer}  [min_{i != j} L_i - L_j]  > 0  (float LP, HiGHS)."""
import itertools, sys, time
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U
E=list(itertools.combinations(range(5),2)); Z=np.array([[(m>>i)&1 for i in range(5)] for m in range(32)])
PAIR=np.array([[Z[m,i]*Z[m,j] for (i,j) in E] for m in range(32)])
C=np.array([[float(x) for x in c.coefficients()] for c in U.make_u248()])
LP=C[:,[0]]+C[:,1:6]@Z.T          # 248 x 32  (pi part, constant absorbed)
LB=C[:,6:]                         # 248 x 10
def base(S):
    A=[];b=[];Aeq=[np.r_[np.ones(32),np.zeros(11)]];beq=[1.0]
    for e,(i,j) in enumerate(E):
        z=np.zeros(10); z[e]=-1
        if e in S: Aeq.append(np.r_[PAIR[:,e],z,0]); beq.append(0.0)
        else: A.append(np.r_[PAIR[:,e],z,0]); b.append(0.0)
        for k in (i,j):
            z=np.zeros(10); z[e]=1; A.append(np.r_[-Z[:,k],z,0]); b.append(0.0)
    for k in range(5):
        A.append(np.r_[Z[:,k],np.zeros(10),0]); b.append(0.5)
        A.append(np.r_[Z.sum(1)-3*Z[:,k],np.zeros(10),0]); b.append(1.0)
    return np.array(A),np.array(b),np.array(Aeq),np.array(beq)
idx={e:k for k,e in enumerate(E)}
layers={'all (Omega_R)':set(),'disjoint M={01,23}':set(range(10))-{idx[(0,1)],idx[(2,3)]},
        'adjacent M={01,02}':set(range(10))-{idx[(0,1)],idx[(0,2)]},'single M={01}':set(range(10))-{idx[(0,1)]},'all_saturated':set(range(10))}
for name,S in layers.items():
    T=time.time(); A,b,Aeq,beq=base(S); nec=0; amb=0
    for j in range(248):
        o=[i for i in range(248) if i!=j]
        R=np.c_[-(LP[o]-LP[j]),-(LB[o]-LB[j]),np.ones(len(o))]
        res=linprog(np.r_[np.zeros(42),-1.0],A_ub=np.vstack([A,R]),b_ub=np.r_[b,np.zeros(len(o))],A_eq=Aeq,b_eq=beq,
                    bounds=[(0,None)]*42+[(None,None)],method='highs')
        t=-res.fun
        if t>1e-7: nec+=1
        elif t>-1e-7: amb+=1
    print(f'{name:22s}: irredundant terms {nec}  (|t|<1e-7 ambiguous: {amb})  [{time.time()-T:.0f}s]',flush=True)
