"""Census: S5-orbits of all vertices of the 60 per-rule dual polyhedra P_F (signed threshold rules, exact cdd enumeration
from d2_98_review/enum_vertices.py) vs. the 248 certificates of U248 and the 11 polyhedra of the report's Lemma 2."""
import json, itertools, sys, time
from fractions import Fraction as F
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U   # report module (read before use)
PERMS=list(itertools.permutations(range(5)))
def canon(v): return min(U.apply_perm(v,p) for p in PERMS)
T=time.time()
u248={canon(c.coefficients()) for c in U.make_u248()}
print('U248 S5-orbits:',len(u248))
data=json.load(open('rule_vertices.json'))
orb={}
for r in data:
    for v in r['vertices']:
        v=tuple(F(x) for x in v); c=canon(v)
        orb.setdefault(c,set()).add(tuple(r['w']))
print('total vertices',sum(len(r['vertices']) for r in data),' distinct S5 orbits',len(orb),' of which in U248:',sum(c in u248 for c in orb),f'[{time.time()-T:.0f}s]')
missing=[c for c in u248 if c not in orb]
print('U248 orbits that are NOT vertices of any P_F:',len(missing))
json.dump([[str(x) for x in c] for c in orb if c not in u248],open('non_u248_orbits.json','w'))
