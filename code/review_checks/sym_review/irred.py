"""Irredundant terms of the reduced family U_det for K = 5..17 (LP over the polygon). Usage: python3 -I irred.py"""
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.argv=[sys.argv[0],'none']
import check_s31 as C
from scipy.optimize import linprog
from fractions import Fraction as Fr
def names(K):
    out={}
    out[(Fr(0),Fr(1),Fr(0))]='dict'; out[(Fr(0),Fr(0),Fr(3))]='3b'
    for s in range(1,(K-1)//2+1):
        for t in range(0,s): out.setdefault(C.odd(s,t),f'B{2*s+1},{t}')
    for s in range(2,K//2+1):
        for t in range(0,s-1):
            if s-1-t>0: out.setdefault(C.even(s,t),f'E{2*s},{t}')
            out.setdefault(C.anch(s,t),f'A{2*s},{t}')
    return out
for K in [5,6,7,8,9,10,11,12,13,14,16,17]:
    T=C.reduced_family(K); nm=names(K); keep=[]
    for i,(a0,B0,C0) in enumerate(T):
        A=[[-1,1,0]]; bu=[0]   # b <= q
        for r in range(K//2+1): A.append([2*r/(K-1),-1,0]); bu.append(r*(r+1)/(K*(K-1)))
        for j,(a,B,Cc) in enumerate(T):
            if j!=i: A.append([float(B0-B),float(C0-Cc),1]); bu.append(float(a-a0))
        r=linprog([0,0,-1],A_ub=A,b_ub=bu,bounds=[(0,.5),(0,.5),(None,1)],method='highs')
        if r.status==0 and -r.fun>1e-9: keep.append(nm.get((a0,B0,C0),'?'))
    print(f'K={K}: reduced family {len(T)} terms, irredundant {len(keep)}: {keep}')
