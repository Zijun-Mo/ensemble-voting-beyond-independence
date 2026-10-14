"""Exact K=5 Lipschitz constants for result 10.
L5  = max over all vertices of the dual polyhedra D_F of the 1882 signed threshold rules (60 S5-classes, exact cdd table
      d2_98_review/rule_vertices.json) of sum|beta| + sum gamma.
L05 = max l1 norm over the 27 zero-intercept atoms of result 3(c) (dictators, triangles, pentagons).
L248= max l1 norm over the 248 certificates of U248 (result 8): constant of the closed-form K=5 certificate C248.
Usage: python3 -I constants_K5.py <rule_vertices.json> <dir with S26 module>"""
import json, sys
from fractions import Fraction as F
d=json.load(open(sys.argv[1])); sys.path.insert(0,sys.argv[2]); import D2_K5_U248_completeness_proof_and_verifier as U
l1=lambda v: sum(abs(x) for x in v[1:6])+sum(v[6:])
best=max((l1([F(x) for x in v]),r['w']) for r in d for v in r['vertices'])
c=[x.coefficients() for x in U.make_u248()]
print('rules:',sum(r['orbit'] for r in d),' L5 =',best[0],'(attained for rule weights',best[1],')')
print('L05 =',max(l1(v) for v in c if v[0]==0),'  L248 =',max(l1(v) for v in c))
