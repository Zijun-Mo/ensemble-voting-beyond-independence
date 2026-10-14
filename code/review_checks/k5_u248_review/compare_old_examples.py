"""Which U248 orbit (if any) are the 3 new certificates found in the 4th review (d2_98_review/examples_1_0.json)?
Usage: python3 -I compare_old_examples.py <dir with S26 module> <path to examples_1_0.json>"""
import sys, json, itertools
from fractions import Fraction as F
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U
P=list(itertools.permutations(range(5))); canon=lambda v: min(U.apply_perm(v,p) for p in P)
lab={}
for c in U.make_u248(): lab.setdefault(canon(c.coefficients()), c.label.split(':')[0])
for d in json.load(open(sys.argv[2])):
    print(lab.get(canon(tuple(F(x) for x in d['cert'])), 'NOT in U248 (dominated, by completeness)'))
