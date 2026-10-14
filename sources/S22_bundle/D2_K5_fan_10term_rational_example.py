#!/usr/bin/env python3
"""Exact rational strict-R fan example: 9-term failure and mixed/strict gap.

Self-contained standard-library check. Deterministic optimality is certified
by the separate 10-term exact-proof program.
"""
from fractions import Fraction as F
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
q=tuple([F(599,2000),F(199,400),F(199,400),F(599,2000),F(599,2000)])
b=tuple([F(797,8000),F(797,8000),F(1,1600),F(1,1600),
         F(797,8000),F(1,1600),F(1,1600)])
qmin=min(q); R=sum(q)-3*qmin
assert R==F(199,200)<1 and all(0<z<F(1,2) for z in q)
assert all(0<=b[k]<min(q[i],q[j]) for k,(i,j) in enumerate(E))
cycles=[(0,1,4),(0,3,5),(1,2,6),(2,3,4,5,6)]
costs=[sum(b[k] for k in cycle) for cycle in cycles]
old=min(*q,*costs)
star=1-q[0]-q[1]-q[2]+sum(b)
assert old==F(807,8000) and star==F(11,1600) and star<old
mixedQ=F(1,2)*(1-q[0]-q[1]-q[2]+b[0]+b[1]+b[4])+b[2]+b[3]+b[5]+b[6]
assert mixedQ==F(3,640)<star

pi_num={1:F(1639,16000),2:F(7,1600),3:F(1549,16000),
        4:F(7,1600),5:F(1549,16000),6:F(1559,16000),
        7:F(7,3200),10:F(4727,16000),11:F(1,1600),
        20:F(4727,16000),21:F(1,1600),24:F(7,3200),
        26:F(1,1600),28:F(1,1600)}
pi=[pi_num.get(m,F(0)) for m in range(32)]
assert sum(pi)==1
assert [sum(pi[m] for m in range(32) if m>>i&1) for i in range(5)]==list(q)
for k,(i,j) in enumerate(E):
    assert sum(pi[m] for m in range(32) if (m>>i&1) and (m>>j&1))==b[k]
nu={3:F(1,1600),5:F(1,1600),7:F(7,3200),
    10:F(1,1600),11:F(1,1600)}
assert sum(nu.values())==F(3,640)
for m,t in nu.items():
    assert 0<=t<=pi[m] and t<=pi[31-m]

print('PASS: strict-R interior fan instance, all moments exact.')
print('Old nine-term candidate =',old,'=',float(old))
print('New deterministic ten-term optimum =',star,'=',float(star))
print('Randomized minimax optimum =',mixedQ,'=',float(mixedQ))
print('Deterministic versus randomized gap =',star-mixedQ)
print('Independent common-adversary primal certificate has symmetric half-mass',sum(nu.values()))
