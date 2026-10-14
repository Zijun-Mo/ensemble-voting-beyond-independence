"""Exact formulation (numerical LP solver) of T_5 for the 7-edge fan graph.

Local-to-global theorem is algebraic: chordal clique tree {012,014,023};
its locally consistent symmetric and residual marginals glue globally.
Only 18 LP variables: t, 7 symmetric pair intersections u_e,
7 residual intersections v_e, and 3 residual triangle intersections z_C.

This script compares the 18-variable LP against the original 48-variable
LP on randomized *feasible* moment caps. Numeric equality tests are
implementation checks, not the mathematical proof of clique gluing.
"""
import numpy as np
from scipy.optimize import linprog
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
TRI=[(0,1,2),(0,1,4),(0,2,3)]
X=np.array([[m>>i&1 for i in range(5)] for m in range(32)],dtype=float)
P=np.array([X[:,i]*X[:,j] for i,j in E]).T

def clique_lp(q,b):
 n=18
 C=[]; d=[]
 def add(coeff,rhs):
  row=np.zeros(n)
  for k,v in coeff.items(): row[k]=v
  C.append(row);d.append(rhs)
 add({0:2},1)
 for i in range(5):add({0:1},q[i])
 for k in range(7):add({1+k:1,8+k:1},b[k])
 for h,tri in enumerate(TRI):
  ids=[E.index(e) for e in __import__('itertools').combinations(tri,2)]
  # symmetric local triangle: sum u_e>=t, 2u_e >= sum u -t
  add({0:1,**{1+k:-1 for k in ids}},0)
  for k in ids:
   row={0:-1,1+k:-1}
   for j in ids:
    row[1+j]=row.get(1+j,0)+1
   row[1+k]-=1
   # equivalent sum(other 2) - u_k -t <= 0
   add(row,0)
  z=15+h
  # residual local: pair-only atoms p110, etc
  for k in ids:add({z:1,8+k:-1},0)
  # residual singleton atoms
  for i in tri:
   k1,k2=[k for k in ids if i in E[k]]
   add({0:1,8+k1:1,8+k2:1,z:-1},q[i])
  # residual p000=1-sum(q_tri)+t+sum(v)-z >=0
  add({0:-1,**{8+k:-1 for k in ids},z:1},1-sum(q[i] for i in tri))
 # Above p000 RHS sign: -t-sum(v)+z <=1-sum(q): yes
 c=np.zeros(n);c[0]=-1
 res=linprog(c,A_ub=C,b_ub=d,bounds=(0,None),method='highs')
 return res

def full_lp(q,b):
 from scipy.optimize import linprog
 Aeq=np.zeros((6,48)); Aeq[0,:32]=1;Aeq[1:,:32]=X.T
 Aub=np.zeros((39,48)); Aub[:7,:32]=P.T
 for j in range(16):
  for m in (j,31-j):Aub[7+m,m]=-1; Aub[7+m,32+j]=1
 c=np.zeros(48);c[32:]=-1
 return linprog(c,A_ub=Aub,b_ub=np.r_[b,np.zeros(32)],A_eq=Aeq,b_eq=np.r_[1,q],bounds=(0,None),method='highs')

if __name__=='__main__':
 rng=np.random.default_rng(17);passed=0
 for trial in range(120):
  p=rng.uniform(.05,.49,size=5)
  pi=np.prod(np.where(X,p,1-p),axis=1)
  q=X.T@pi;actual=P.T@pi
  b=actual+(np.minimum(q[[i for i,j in E]],q[[j for i,j in E]])-actual)*rng.uniform(0,1,size=7)
  x=clique_lp(q,b);y=full_lp(q,b)
  assert x.success and y.success,(x.message,y.message)
  assert abs(x.fun-y.fun)<1e-8,(x.fun,y.fun)
  passed+=1
 print('PASSED:',passed,'feasible instances, reduced 18-variable and full 48-variable LP agree.')
