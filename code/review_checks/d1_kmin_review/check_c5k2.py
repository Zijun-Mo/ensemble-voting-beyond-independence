import itertools, sympy as sp, networkx as nx
from fractions import Fraction as F
import sys
sys.path.insert(0,'.')
def odd_cycles(n):
    # all simple odd cycles of K_n as edge sets
    cyc=[]
    V=range(n)
    for L in range(3,n+1,2):
        for S in itertools.combinations(V,L):
            s0=S[0]
            for perm in itertools.permutations(S[1:]):
                if perm[0]<perm[-1]:
                    seq=(s0,)+perm
                    cyc.append([tuple(sorted((seq[i],seq[(i+1)%L]))) for i in range(L)])
    return cyc
def check(n,u,name):
    E=list(itertools.combinations(range(n),2)); idx={e:i for i,e in enumerate(E)}
    cyc=odd_cycles(n)
    mins=min(sum(u[e] for e in c) for c in cyc)
    tight=[c for c in cyc if sum(u[e] for e in c)==1]
    rows=[]
    for c in tight:
        r=[0]*len(E)
        for e in c: r[idx[e]]+=1
        rows.append(r)
    for e in E:
        if u[e]==0:
            r=[0]*len(E); r[idx[e]]=1; rows.append(r)
    rank=sp.Matrix(rows).rank() if rows else 0
    print(name,'n',n,'min odd-cycle sum',mins,'tight cycles',len(tight),'rank',rank,'/',len(E))
    return E
def c5k2():
    u={}
    for e in itertools.combinations(range(7),2):
        a,b=e
        if a<5 and b<5:
            u[e]=F(1,5) if (b-a)%5 in (1,4) else F(3,5)
        elif a<5: u[e]=F(2,5)
        else: u[e]=F(1,5)
    return u
u7=c5k2(); check(7,u7,'C5+K2')
def lift(u,w):
    v={}
    for e in itertools.combinations(range(8),2):
        a,b=e
        if b<7: v[e]=u[e]
        elif a==w: v[e]=F(0)
        else: v[e]=1-u[tuple(sorted((a,w)))]
    return v
uA=lift(u7,6); check(8,uA,'LiftA')
uB=lift(u7,0); check(8,uB,'LiftB')
import pickle; pickle.dump((u7,uA,uB),open('c5k2.pkl','wb'))
