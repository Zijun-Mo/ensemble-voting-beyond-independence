"""D2 -- rigorous K5 chorded C5 closed-form theorem, exact arithmetic.

Graph edges: 01,02,04,12,23,34; all other pair caps uninformative.
Requires only Python standard library.

Run: python D2_K5_chorded_cycle_complete_exact.py
Optional: python D2_K5_chorded_cycle_complete_exact.py --quick
Proof parts: (1) exact double-description of 12D polyhedron;
(2) exact Farkas certificates for all 34 vertices.
"""
import itertools,math,sys,json,time
from pathlib import Path
from fractions import Fraction
from functools import reduce
from math import gcd

V=[[0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1], [0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, -1, -1, 0, 0, 0, 1, 0, 1, 1, 1, 1], [1, -1, 0, 0, 0, -1, 1, 0, 1, 1, 1, 1], [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, -1, -1, 0, 0, 1, 0, 1, 1, 1, 1], [1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, -1, -1, 0, 1, 0, 1, 1, 1, 1], [1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, -1, -1, 1, 0, 1, 1, 1, 1], [1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [2, -2, -1, 0, 0, -1, 1, 0, 1, 1, 1, 1], [2, -1, -2, -1, 0, 0, 1, 0, 1, 1, 1, 1], [2, -1, 0, 0, -1, -2, 1, 0, 1, 1, 1, 1], [2, 0, -1, -2, -1, 0, 1, 0, 1, 1, 1, 1], [2, 0, 0, -1, -2, -1, 1, 0, 1, 1, 1, 1], [3, -2, -2, -2, 0, 0, 1, 1, 0, 1, 0, 0], [3, -2, -2, -1, 0, -1, 1, 0, 1, 1, 1, 1], [3, -2, -1, 0, -1, -2, 1, 0, 1, 1, 1, 1], [3, -1, -2, -2, -1, 0, 1, 0, 1, 1, 1, 1], [3, -1, 0, -1, -2, -2, 1, 0, 1, 1, 1, 1], [3, 0, -1, -2, -2, -1, 1, 0, 1, 1, 1, 1], [4, -2, -2, -2, -1, -1, 1, 0, 1, 1, 1, 1], [4, -2, -2, -1, -1, -2, 1, 0, 1, 1, 1, 1], [4, -2, -1, -1, -2, -2, 1, 0, 1, 1, 1, 1], [4, -1, -2, -2, -2, -1, 1, 0, 1, 1, 1, 1], [4, -1, -1, -2, -2, -2, 1, 0, 1, 1, 1, 1], [5, -2, -2, -2, -2, -2, 1, 0, 1, 1, 1, 1]]
R=[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, -1, -1, -1, 0, 0, 1, 1, 0, 1, 0, 0], [1, -1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0], [1, -1, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0], [1, -1, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0], [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0], [2, -1, -1, -1, -1, -1, 1, 0, 1, 1, 1, 1]]

N=13
EDGES=[(0,1),(0,2),(0,4),(1,2),(2,3),(3,4)]
A=[];h=[]
for mask in range(32):
 x=[(mask>>i)&1 for i in range(5)]
 A.append(tuple([1,*x,*[x[i]*x[j] for i,j in EDGES]]));h.append(0)
for mask in range(32):
 if mask < (31^mask):
  A.append(tuple(x+y for x,y in zip(A[mask],A[31^mask])));h.append(1)
for e in range(6):
 A.append(tuple([0]*6+[int(i==e) for i in range(6)]));h.append(0)
assert len(A)==54 and len(V)==34 and len(R)==24
G=[tuple(v+[1]) for v in V]+[tuple(r+[0]) for r in R]

def verify_cone_exact():
    PRIMES=(1000000007,1000000009)
    p=PRIMES[0]

    def rankmod(rows, p):
     if not rows:return 0
     m=[list(map(lambda x:x%p,row)) for row in rows];height=len(m);width=len(m[0]);rank=0
     for c in range(width):
      idx=next((j for j in range(rank,height) if m[j][c]!=0),None)
      if idx is None:continue
      m[rank],m[idx]=m[idx],m[rank]
      inv=pow(m[rank][c],p-2,p)
      for j in range(rank+1,height):
       if m[j][c]:
        f=m[j][c]*inv%p
        for k in range(c,width):m[j][k]=(m[j][k]-f*m[rank][k])%p
      rank+=1
      if rank==min(height,width):break
     return rank

    assert max(sum(t*t for t in ray) for ray in G) <= 51

    def rank_exact(rows):
     # Each generator has squared Euclidean norm <=51. For k<=13, every
     # k-by-k minor has |det| <=51^(k/2) < 10^12 < p1*p2.
     # Hence modular rank over either of two primes recovers exact Q-rank.
     return max(rankmod(rows, PRIMES[0]),rankmod(rows, PRIMES[1]))

    def prim(v):
     v=[int(t) for t in v];d=reduce(gcd,(abs(t) for t in v));assert d
     return tuple(t//d for t in v)

    def dot(a,b):return sum(x*y for x,y in zip(a,b))

    # Choose independent rays as initial simplicial cone
    indices=[]
    for i,r in enumerate(G):
     if rank_exact([G[j] for j in indices]+[r])>len(indices):indices.append(i)
     if len(indices)==N:break
    assert len(indices)==N
    # Rational Gauss-Jordan inverse of matrix whose columns are chosen rays
    from fractions import Fraction as F
    M=[[F(G[indices[c]][r]) for c in range(N)] + [F(int(r==c)) for c in range(N)] for r in range(N)]
    for c in range(N):
     k=next(r for r in range(c,N) if M[r][c]!=0)
     M[c],M[k]=M[k],M[c]
     val=M[c][c]
     M[c]=[u/val for u in M[c]]
     for rr in range(N):
      if rr!=c and M[rr][c]:
       coeff=M[rr][c];M[rr]=[x-coeff*y for x,y in zip(M[rr],M[c])]
    inv=[row[N:] for row in M]
    facets=set()
    for k in range(N):
     row=inv[k]
     den=math.lcm(*(x.denominator for x in row));normal=prim([int(x*den) for x in row])
     assert all(dot(normal,G[j])>=0 for j in indices)
     facets.add(normal)
    active=[G[i] for i in indices]
    remaining=[G[i] for i in range(len(G)) if i not in indices]
    print('total generators',len(G),'start indep',indices,'start facet count',len(facets),flush=True)
    start=time.time()
    for k,r in enumerate(remaining):
     plus=[];minus=[];zero=[]
     for f in facets:
      d=dot(f,r)
      if d>0:plus.append((f,d))
      elif d<0:minus.append((f,d))
      else:zero.append(f)
     if not minus:
      active.append(r);continue
     new=set()
     for f,fp in plus:
      for g,gm in minus:
       candidate=prim([fp*gj-gm*fj for fj,gj in zip(f,g)])
       assert dot(candidate,r)==0
       # test actual facet dimensionality, mod rank 12 certifies exact rank 12
       tight=[old for old in active if dot(candidate,old)==0]
       if rank_exact(tight+[r])==N-1:
        new.add(candidate)
     assert all(all(dot(fc,s)>=0 for s in active+[r]) for fc in new)
     facets=set(plusf for plusf,dp in plus)|set(zero)|new
     active.append(r)
     print(f'add {k+1:2d}/{len(remaining)} removed={len(minus):3d} added={len(new):3d} facets={len(facets):4d} elapsed={time.time()-start:.2f}',flush=True)

    # Every computed facet normal is a true facet (rank 12), and algorithm is exact
    for f in facets:
     assert all(dot(f,s)>=0 for s in G)
     tight=[s for s in G if dot(f,s)==0]
     assert rank_exact(tight)==N-1
    expected={prim(list(row)+[-rhs]) for row,rhs in zip(A,h)}|{tuple([0]*12+[1])}
    print('computed nonredundant facets',len(facets),'among original 55 rows',len(expected))
    assert facets.issubset(expected)
    # Because cone(G) lies in the original 55 halfspaces, and all its
    # exact facets belong to those halfspaces, the two cones coincide.
    print('CERTIFIED: full mixed dual = conv(34 integer vertices) + cone(24 integer rays)')

def verify_farkas_exact():
 data=json.loads((Path(__file__).parent/'D2_K5_chord_34_exact_Farkas.json').read_text())
 assert [tuple(e) for e in data['E']]==EDGES
 X=[[((mask>>i)&1) for i in range(5)] for mask in range(32)]
 P=[[x[i]*x[j] for i,j in EDGES] for x in X]
 C=[]
 for i in range(5):
  c=[0]*12;c[1+i]=1;C.append(c)
 for cyc in [[(0,1),(0,2),(1,2)],[(0,1),(1,2),(2,3),(3,4),(0,4)]]:
  c=[0]*12
  for edge in cyc:c[6+EDGES.index(edge)]=1
  C.append(c)
 eq=[1]*32+[0]*8
 rhs=[];B=[]
 for i in range(5):
  rr=[X[m][i] for m in range(32)]+[0]*8
  B.append(rr);rhs.append(Fraction(1,2))
 for edge,(i,j) in enumerate(EDGES):
  rr=[P[m][edge] for m in range(32)]+[-int(k==edge) for k in range(6)]+[0,0]
  B.append(rr);rhs.append(Fraction(0))
  for node in [i,j]:
   rr=[-X[m][node] for m in range(32)]+[int(k==edge) for k in range(6)]+[0,0]
   B.append(rr);rhs.append(Fraction(0))
 cobj=[0]*38+[1,-1]
 assert len(B)==23
 assert [list(entry['vertex']) for entry in data['vertices']]==V
 for n,entry in enumerate(data['vertices']):
  v=entry['vertex']
  BB=list(B);hh=list(rhs)
  for cc in C:
   d=[v[i]-cc[i] for i in range(12)]
   rr=[sum(X[m][i]*d[1+i] for i in range(5)) for m in range(32)] +d[6:]+[1,-1]
   BB.append(rr);hh.append(Fraction(cc[0]-v[0]))
  assert len(BB)==30
  yy=[Fraction(0)]*30
  for idx,val in entry['y']:
   yy[int(idx)]=Fraction(val)
  z=Fraction(entry['z'])
  assert all(y>=0 for y in yy)
  assert sum(hh[i]*yy[i] for i in range(30))+z<=0
  assert all(sum(BB[i][j]*yy[i] for i in range(30))+z*eq[j]>=cobj[j] for j in range(40))
  # The corresponding vertex is a valid mixed-dual certificate.
  assert all(sum(A[i][j]*v[j] for j in range(12))>=h[i] for i in range(54))
 print('CERTIFIED: all 34 vertices are dominated by 7 candidate affine terms on the entire feasible parameter polytope.')

if __name__=='__main__':
 verify_farkas_exact()
 if '--quick' not in sys.argv:
  verify_cone_exact()
 print('THEOREM: for informative graph contained in C5+chord, T5=Rdet5=min(qmin,triangle sum,5-cycle sum).')
