"""Exact rational checks for the K=4 structural classification and its K=5 lift.
No random sampling, optimization solver, or floating point arithmetic.
Requires Python >= 3.10. Run as a standalone file.
"""
from fractions import Fraction as F
from itertools import combinations, permutations

N4=4
states4=[tuple((m>>i)&1 for i in range(4)) for m in range(16)]
compl4=lambda m: m ^ 15

def truth_of_w(w):
    assert all(sum(v*(2*x[i]-1) for i,v in enumerate(w)) != 0 for x in states4)
    return tuple(int(sum(v*(2*x[i]-1) for i,v in enumerate(w)) > 0) for x in states4)

# The 12 monotone self-dual Boolean functions: dictators, 3-core majority,
# and four-anchor (2,1,1,1) majority. All truth tables are distinct.
mono={}
for i in range(4):
    w=[0]*4;w[i]=1
    mono[truth_of_w(w)]=('dict',i)
for i in range(4):
    w=[1]*4;w[i]=0
    mono[truth_of_w(w)]=('majority3',i)
for i in range(4):
    w=[1]*4;w[i]=2
    mono[truth_of_w(w)]=('anchor4',i)
assert len(mono)==12

# Exhaust all 2^8=256 self-dual rules and *independently* check that
# precisely 12 are monotone, confirming the complete classification.
reps=[m for m in range(16) if m<compl4(m)]
monotone_tables=set()
for mask in range(256):
    f=[0]*16
    for j,m in enumerate(reps):
        f[m]=(mask>>j)&1
        f[compl4(m)]=1-f[m]
    is_monotone=all(f[m]<=f[n] for m in range(16)
                    for n in range(16) if (m&n)==m)
    if is_monotone:monotone_tables.add(tuple(f))
assert len(monotone_tables)==12 and monotone_tables==set(mono)

# Input coordinate sign flips of the 12 basic weighted rules generate
# *all* linear threshold rules (proved by the sign-normalization lemma).
signed={}
for f,(kind,i) in mono.items():
    for flipmask in range(16):
        g=tuple(f[m^flipmask] for m in range(16))
        signed[g]=kind
assert len(signed)==104
assert {t:sum(k==t for k in signed.values()) for t in ('dict','majority3','anchor4')}=={
    'dict':8,'majority3':32,'anchor4':64}

# Construct a genuinely positive-intercept 4-model certificate for anchor=0.
# Q=1-sum_i x_i + sum_{j<k; j,k>0}x_j*x_k + 2*sum_{j>0}x_0*x_j.
for x in states4:
    Q=1-sum(x)+sum(x[j]*x[k] for j,k in combinations((1,2,3),2)) \
      +2*sum(x[0]*x[j] for j in (1,2,3))
    fail=int(2*x[0]+sum(x[1:])>=3)
    assert Q>=fail and Q>=0

# Exact one-parameter zero-risk family, p=2/5 given for a finite check.
p=F(2,5)
support4={1:F(2,5)}           # X0=1 and the other 3 models correct
for i,j in combinations((1,2,3),2):
    support4[(1<<i)|(1<<j)]=(1-p)/3
assert sum(support4.values())==1
q4=tuple(sum(v*((mask>>i)&1) for mask,v in support4.items()) for i in range(4))
E4=list(combinations(range(4),2))
b4={e:sum(v*((mask>>e[0])&1)*((mask>>e[1])&1)
         for mask,v in support4.items()) for e in E4}
assert q4==(p,2*(1-p)/3,2*(1-p)/3,2*(1-p)/3)
assert all(b4[(0,j)]==0 for j in (1,2,3))
assert all(b4[(i,j)]==(1-p)/3 for i,j in combinations((1,2,3),2))
B_anchor=1-sum(q4)+sum(b4[e] for e in combinations((1,2,3),2)) \
  +2*sum(b4[(0,j)] for j in (1,2,3))
assert B_anchor==0

# Lift to K=5 by duplicating X0 as X4; this is a feasible joint law in R.
support5={m|(((m>>0)&1)<<4):v for m,v in support4.items()}
q5=[sum(v*((m>>i)&1) for m,v in support5.items()) for i in range(5)]
E5=list(combinations(range(5),2))
b5={e:sum(v*((m>>e[0])&1)*((m>>e[1])&1)
          for m,v in support5.items()) for e in E5}
assert q5==[p,2*(1-p)/3,2*(1-p)/3,2*(1-p)/3,p]
assert sum(q5)-3*min(q5)<1 and all(0<q<F(1,2) for q in q5)
assert all(2*((m>>0)&1)+sum((m>>i)&1 for i in (1,2,3)) < 3
           for m in support5)
assert min(q5)>0 and min(b5.values())==0
# The 10 triangles and 12 Hamilton cycles all have strictly positive cost.
odd_cycles=set()
for T in combinations(range(5),3):
    odd_cycles.add(tuple(sorted(combinations(T,2))))
for perm in permutations(range(1,5)):
    seq=(0,*perm)
    odd_cycles.add(tuple(sorted(tuple(sorted((seq[i],seq[(i+1)%5]))) for i in range(5))))
assert len(odd_cycles)==22
assert min(sum(b5[e] for e in C) for C in odd_cycles)>0
print('PASS K4: 12 monotone self-dual, 104 signed threshold truth tables')
print('PASS Q_anchor(x)>=F_anchor(x), exact on all 16 states')
print('PASS K4: positive-intercept 4-core attains robust failure probability zero')
print('PASS K5: strictly inside R, q_min>0, all odd cycles positive, 4-core risk zero')
print('q5 =',q5, '; worst bound =',B_anchor)
