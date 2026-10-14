import sys, numpy as np, itertools
from fractions import Fraction as F
import os; sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','..','sources','S21_bundle'))
from D2_K5_fan_8var_exact_LP import fan_8var_lp, fan_full_lp
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
q=[F(599,2000),F(995,2000),F(995,2000),F(599,2000),F(599,2000)]; b=[F(x,8000) for x in (797,797,5,5,797,5,5)]
# Q_mix validity over all 32 states, exactly
def Q(x): return F(1,2)*(1-x[0]-x[1]-x[2]+x[0]*x[1]+x[0]*x[2]+x[1]*x[2])+x[0]*x[3]+x[0]*x[4]+x[1]*x[4]+x[2]*x[3]
X=[[m>>i&1 for i in range(5)] for m in range(32)]
print('Q_mix >= 0 everywhere:',all(Q(x)>=0 for x in X),'; Q(x)+Q(1-x) >= 1 everywhere:',all(Q(x)+Q([1-t for t in x])>=1 for x in X))
L=F(1,2)*(1-q[0]-q[1]-q[2]+b[0]+b[1]+b[4])+b[2]+b[3]+b[5]+b[6]
cyc=[b[0]+b[1]+b[4],b[0]+b[3]+b[5],b[1]+b[2]+b[6],b[2]+b[3]+b[4]+b[5]+b[6]]
t0=min(*q,*cyc); Bfan=1-q[0]-q[1]-q[2]+sum(b)
print('cycle costs (x8000):',[c*8000 for c in cyc],' t0 =',t0,' B_fan =',Bfan,' L(Q_mix) =',L)
a=fan_8var_lp(np.array(q,float),np.array(b,float)); f=fan_full_lp(np.array(q,float),np.array(b,float))
print(f'8-variable LP T_H = {-a.fun:.7f}; 48-variable LP T_H = {-f.fun:.7f}; 3/640 = {3/640:.7f}')
