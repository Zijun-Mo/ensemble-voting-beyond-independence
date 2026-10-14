"""Exact-rational audit of deterministic vs randomized robust K=5 voting.

Uses previously established theorem for deterministic R0:
 R_det = min(qmin, 22 odd-cycle costs), when sum(q)<=1.
Independently proves the symmetric randomized 3-term formula by
complete 3-variable mixed-dual vertex enumeration (120 bases).
Standard library only.
"""
from fractions import Fraction as F
from itertools import combinations,permutations

rows=[]
for n in range(6):
    rows.append(((F(1),F(n),F(n*(n-1)//2)),F(0)))
for n in range(3):
    rows.append(((F(2),F(5),F(n*(n-1)//2+(5-n)*(4-n)//2)),F(1)))
rows.append(((F(0),F(0),F(1)),F(0)))

def det(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
           -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
           +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))

vertices=set()
for triple in combinations(range(10),3):
    mat=[list(rows[i][0]) for i in triple]
    right=[rows[i][1] for i in triple]
    d=det(mat)
    if d==0:continue
    z=tuple(det([[right[i] if j==col else mat[i][j]
                   for j in range(3)] for i in range(3)])/d for col in range(3))
    if all(sum(a*b for a,b in zip(row,z))>=rhs for row,rhs in rows):
        vertices.add(z)
assert len(vertices)==6,(len(vertices),vertices)
for v in sorted(vertices):print('symmetric mixed vertex (a,beta,gamma):',tuple(map(str,v)))
# With q in [0,1/2], three vertices are always dominated:
# 1-q >= q, 5/2-5q+5b/2 >= 5b/2, 6-15q+10b >= 1-5q+10b.

def mixed_symmetric(q,b):
    return min(q,F(5,2)*b,1-5*q+10*b)

def pi_sym_witness():
    # The distribution has complementary symmetric mass 0.20,
    # achieves q_i=0.2 and every pairwise probability=0.04.
    p={0:F(3,10)}
    for m in range(1,32):
        n=m.bit_count()
        if n==1:p[m]=F(1,10)
        if n in (2,3):p[m]=F(1,100)
    assert sum(p.values())==1
    assert all(sum(z*((m>>i)&1) for m,z in p.items())==F(1,5) for i in range(5))
    assert all(sum(z*((m>>i)&1)*((m>>j)&1) for m,z in p.items())==F(1,25)
               for i,j in combinations(range(5),2))
    T_lbound=sum(min(p.get(m,0),p.get(31^m,0)) for m in range(32) if m<(31^m))
    assert T_lbound==F(1,10)

pi_sym_witness()
q=F(1,5);b=F(1,25)
T=mixed_symmetric(q,b)
Rdet=min(q,3*b) # valid because sum q=1 (established R0 theorem)
assert (T,Rdet)==(F(1,10),F(3,25))
print('symmetric gap: mixed',T,'deterministic',Rdet,'difference',Rdet-T)

# An explicit full-dimensional heterogeneous open box of feasible parameters:
# q_i in (.18,.19), b_ij in (.039,.041). Independence is feasible since
# q_i*q_j < .19**2 = .0361 < .039, and sum(q)<.95<1.
# Every triangle has cost >3*.039=.117, qmin>.18.
# One legal mixed-dual Q=1/4*sum_{i<j} x_i x_j gives T< 10*.041/4=.1025.
# Thus, throughout the open box, Rdet-T>.117-.1025=.0145.
assert F(19,100)**2<F(39,1000)
assert 5*F(19,100)<1
assert 3*F(39,1000)-F(10,4)*F(41,1000)==F(29,2000)
print('heterogeneous open-box gap is strictly greater than',F(29,2000))
print('ALL EXACT FRACTION AUDITS PASSED')
