#!/usr/bin/env python3
"""Certified fan-graph K=5 deterministic minimax 10-term formula on R.

Run: python D2_K5_fan_10term_exact_verify.py
Requires: Python >=3.10 standard library only.

Input certificates (same folder):
  D2_K5_fan_10term_generators.json
  D2_K5_fan_10term_Farkas_witnesses.json

Verification is EXACT: integer double-description, modular rank certified by
Hadamard's bound, and rational Farkas identities for every one-pair vertex.
"""
import json
from fractions import Fraction as F
from functools import reduce
from math import gcd
from pathlib import Path

BASE=Path(__file__).resolve().parent
EDGES=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
CYCLES=[[(0,1),(0,2),(1,2)],[(0,1),(0,4),(1,4)],
        [(0,2),(0,3),(2,3)],[(0,3),(0,4),(1,2),(1,4),(2,3)]]
X=[[m>>i&1 for i in range(5)] for m in range(32)]
M=[tuple([1,*x,*[x[i]*x[j] for i,j in EDGES]]) for x in X]
D=13
A=[*M]
h=[0]*32
for m in range(16):
    A.append(tuple(M[m][i]+M[31-m][i] for i in range(D)));h.append(1)
for k in range(7):
    A.append(tuple([0]*6+[int(j==k) for j in range(7)]));h.append(0)
assert len(A)==55 and len(h)==55
H=[tuple([*row,-rhs]) for row,rhs in zip(A,h)]+[tuple([0]*13+[1])]
N=14
atlas=json.loads((BASE/'D2_K5_fan_10term_generators.json').read_text())
V=[tuple(map(int,v)) for v in atlas['V']]
R=[tuple(map(int,r)) for r in atlas['R']]
G=sorted(set(V+R))
assert len(V)==43 and len(R)==28 and len(G)==71
assert all(z[-1]>0 for z in V) and all(z[-1]==0 for z in R)
assert all(all(sum(a*b for a,b in zip(row,g))>=0 for row in H) for g in G)


def primitive(v):
    div=reduce(gcd,(abs(int(t)) for t in v))
    assert div>0
    return tuple(int(t)//div for t in v)


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def rankmod(rows,p):
    if not rows:return 0
    mat=[[int(x)%p for x in row] for row in rows]
    nr=len(mat);nc=len(mat[0]);rank=0
    for c in range(nc):
        k=next((i for i in range(rank,nr) if mat[i][c]),None)
        if k is None:continue
        mat[rank],mat[k]=mat[k],mat[rank]
        inv=pow(mat[rank][c],p-2,p)
        for j in range(rank+1,nr):
            if mat[j][c]:
                f=mat[j][c]*inv%p
                for z in range(c,nc):mat[j][z]=(mat[j][z]-f*mat[rank][z])%p
        rank+=1
        if rank==min(nr,nc):break
    return rank

P1=1000000007;P2=1000000009
MAXNORM=max(dot(g,g) for g in G)
assert MAXNORM==179 and MAXNORM**7<P1*P2
# Hadamard: every square minor of <=14 rows has absolute determinant
# <= (max row squared Euclidean norm)^(14/2) < P1*P2.
# Hence rank over Q is max(rank mod P1,rank mod P2).
def rank_exact(rows):return max(rankmod(rows,P1),rankmod(rows,P2))


def verify_full_cone():
    """Exact incremental double-description on the 71 homogeneous rays."""
    indices=[]
    for i,r in enumerate(G):
        if rank_exact([G[j] for j in indices]+[r])>len(indices):indices.append(i)
        if len(indices)==N:break
    assert len(indices)==N
    # exact inverse of independent ray matrix, columns are rays
    mat=[[F(G[indices[c]][r]) for c in range(N)]+[F(int(r==c)) for c in range(N)] for r in range(N)]
    for c in range(N):
        k=next(r for r in range(c,N) if mat[r][c])
        mat[c],mat[k]=mat[k],mat[c]
        pv=mat[c][c];mat[c]=[x/pv for x in mat[c]]
        for j in range(N):
            if j!=c and mat[j][c]:
                val=mat[j][c];mat[j]=[x-val*y for x,y in zip(mat[j],mat[c])]
    import math
    facets=set()
    for k in range(N):
        row=mat[k][N:]
        den=math.lcm(*(z.denominator for z in row))
        ff=primitive([int(z*den) for z in row])
        assert all(dot(ff,G[j])>=0 for j in indices)
        facets.add(ff)
    active=[G[i] for i in indices]
    remaining=[G[i] for i in range(len(G)) if i not in indices]
    for r in remaining:
        pos=[];neg=[];zero=[]
        for f in facets:
            d=dot(f,r)
            if d>0:pos.append((f,d))
            elif d<0:neg.append((f,d))
            else:zero.append(f)
        if not neg:
            active.append(r)
            continue
        new=set()
        for f,fp in pos:
            for g,gm in neg:
                n=primitive([fp*gj-gm*fj for fj,gj in zip(f,g)])
                assert dot(n,r)==0
                tight=[old for old in active if dot(n,old)==0]
                if rank_exact(tight+[r])==N-1:new.add(n)
        assert all(all(dot(fc,s)>=0 for s in active+[r]) for fc in new)
        facets={f for f,_ in pos}|set(zero)|new
        active.append(r)
    assert len(facets)==49
    expected={primitive(row) for row in H}
    assert facets<=expected,('unaccounted facets',facets-expected)
    print('PASS cone: 71 generators, 49 exact facets, all from the original 56 inequalities')


def poly_at(Q,mask):return sum(F(z)*Q[j] for j,z in enumerate(M[mask]))


def candidates(branch):
    """Finite superset of vertices after one pair's disjunctive cut."""
    vv=[tuple(F(z,v[-1]) for z in v[:-1]) for v in V]
    rr=[tuple(F(z) for z in r[:-1]) for r in R]
    val=[poly_at(v,branch) for v in vv]
    good=[i for i,x in enumerate(val) if x>=1]
    bad=[i for i,x in enumerate(val) if x<1]
    greater=[i for i,x in enumerate(val) if x>1]
    out={vv[i] for i in good}
    for i in bad:
        for j in greater:
            theta=(1-val[i])/(val[j]-val[i])
            assert 0<theta<1
            out.add(tuple((1-theta)*vv[i][k]+theta*vv[j][k] for k in range(D)))
        for r in rr:
            slope=poly_at(r,branch)
            if slope>0:
                t=(1-val[i])/slope
                out.add(tuple(vv[i][k]+t*r[k] for k in range(D)))
    assert len(out)==239
    for Q in out:
        assert poly_at(Q,branch)>=1
        assert all(poly_at(Q,m)>=0 for m in range(32))
        assert all(poly_at(Q,m)+poly_at(Q,31-m)>=1 for m in range(16))
        assert all(z>=0 for z in Q[6:])
    return out

STAR=tuple(map(F,[1,-1,-1,-1,0,0,*([1]*7)]))
CERT=[]
for i in range(5):
    z=[0]*D;z[1+i]=1;CERT.append(tuple(map(F,z)))
for cyc in CYCLES:
    z=[0]*D
    for edge in cyc:z[6+EDGES.index(edge)]=1
    CERT.append(tuple(map(F,z)))
CERT.append(STAR)
assert len(CERT)==10


def verify_upper_bounds():
    for Q in CERT:
        assert all(poly_at(Q,m)>=0 for m in range(32))
        assert all(max(poly_at(Q,m),poly_at(Q,31-m))>=1 for m in range(16))
        assert all(z>=0 for z in Q[6:])
    print('PASS upper: 10 quadratic certificates orient to deterministic self-dual rules')


def parameter_rows():
    """R's feasible moment parameter polytope, plus t<=10 candidates."""
    B=[];rhs=[]
    def put(pi=None,b=None,t=0,val=0):
        row=[F(0)]*40
        if pi is not None:
            for m,z in enumerate(pi):row[m]=F(z)
        if b is not None:
            for k,z in enumerate(b):row[32+k]=F(z)
        row[39]=F(t);B.append(tuple(row));rhs.append(F(val))
    for k,(i,j) in enumerate(EDGES):
        b=[0]*7;b[k]=-1
        put(pi=[X[m][i]*X[m][j] for m in range(32)],b=b)
        for v in (i,j):
            b=[0]*7;b[k]=1
            put(pi=[-X[m][v] for m in range(32)],b=b)
    for i in range(5):put(pi=[X[m][i] for m in range(32)],val=F(1,2))
    for Q in CERT:
        put(pi=[-sum(Q[1+i]*X[m][i] for i in range(5)) for m in range(32)],
            b=[-v for v in Q[6:]],t=1,val=Q[0])
    for j in range(5):
        put(pi=[sum(X[m])-3*X[m][j] for m in range(32)],val=1)
    assert len(B)==41
    return B,rhs


def verify_farkas_exact():
    B,hb=parameter_rows()
    entries=json.loads((BASE/'D2_K5_fan_10term_Farkas_witnesses.json').read_text())['entries']
    assert len(entries)==478
    for branch in (7,24):
        listed=[z for z in entries if z['branch']==branch]
        cc=candidates(branch)
        assert len(listed)==239
        assert {tuple(F(s) for s in z['Q']) for z in listed}==cc
    for ent in entries:
        Q=tuple(F(s) for s in ent['Q'])
        yy=[F(0)]*len(B)
        for i,value in ent['y']:
            i=int(i);assert yy[i]==0
            yy[i]=F(value)
        lam=F(ent['lambda'])
        assert min(yy)>=0
        c=[sum(Q[1+i]*X[m][i] for i in range(5)) for m in range(32)]+list(Q[6:])+[F(-1)]
        for j in range(40):
            z=c[j]+sum(yy[i]*B[i][j] for i in range(len(B)))
            if j<32:z-=lam
            assert z>=0,('stationarity',ent['branch'],j,z)
        assert Q[0]+lam-sum(yy[i]*hb[i] for i in range(len(B)))>=0
    print('PASS Farkas: 478 rational witnesses (239 per orientation) over entire feasible R')

if __name__=='__main__':
    verify_upper_bounds()
    verify_full_cone()
    verify_farkas_exact()
    print('THEOREM: for any feasible q in [0,1/2]^5 with sum(q)-3*min(q)<=1,')
    print('R_det,5 = min(q_i, three fan triangles, fan C5, 1-q0-q1-q2+sum(fan b)).')
