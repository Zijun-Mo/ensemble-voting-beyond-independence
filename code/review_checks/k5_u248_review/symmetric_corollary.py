"""Result 8(c): the symmetric K=5 case (q_i=q, b_ij=b) lies entirely in closure(R) (sum q - 3 q_min = 2q <= 1),
so U248 must reduce to min{q, 3b, 1-5q+10b}. Exact check on random rational feasible symmetric points + float sweep.
Usage: python3 -I symmetric_corollary.py <dir with S26 module> <dir with exact_core.py>"""
import sys, numpy as np
from fractions import Fraction as F
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[2])
import D2_K5_U248_completeness_proof_and_verifier as U
from exact_core import feasible
certs=U.make_u248(); rng=np.random.default_rng(1); n=0; bad=0
while n<3000:
    q=F(int(rng.integers(1,500)),1000); b=F(int(rng.integers(0,1000)),1000)*q
    if not feasible(np.full(5,float(q)),np.full(10,float(b))): continue
    n+=1; u=min(c.bound([q]*5,[b]*10) for c in certs)
    bad+= u!=min(q,3*b,1-5*q+10*b)
print(n,'feasible symmetric rational points; exact mismatches between U248 and min{q,3b,1-5q+10b}:',bad)
