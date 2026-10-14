"""Integer-only check: 1,882 signed threshold rules and 60 S5-orbits."""
from itertools import permutations, product
CANON=[(1,0,0,0,0),(3,1,1,1,1),(2,1,1,1,0),
       (3,2,2,1,1),(1,1,1,0,0),(2,2,1,1,1),(1,1,1,1,1)]
REP=[s for s in range(32) if not(s&1)]
PERMS=list(permutations(range(5)))

def bits_from_w(w):
    f=[]
    for s in range(32):
        dot=sum(w[i]*(1 if s&(1<<i) else -1) for i in range(5))
        assert dot!=0
        f.append(int(dot>0))
    assert all(f[s]+f[s^31]==1 for s in range(32))
    return sum(f[s]<<j for j,s in enumerate(REP))

def permute_code(code,perm):
    f=[0]*32
    for j,s in enumerate(REP):
        f[s]=(code>>j)&1
        f[s^31]=1-f[s]
    out=[0]*32
    for s in range(32):
        image=sum(((s>>i)&1)<<perm[i] for i in range(5))
        out[image]=f[s]
    return sum(out[s]<<j for j,s in enumerate(REP))

rules=set()
for w0 in CANON:
    for w in set(permutations(w0)):
        for sg in product((-1,1), repeat=5):
            rules.add(bits_from_w(tuple(w[i]*sg[i] for i in range(5))))
assert len(rules)==1882
seen=set(); orbits=0
for code in rules:
    if code in seen: continue
    orb={permute_code(code,p) for p in PERMS}
    assert orb.issubset(rules)
    seen|=orb
    orbits+=1
assert seen==rules and orbits==60
print('signed self-dual threshold rules:',len(rules),'S5 orbits:',orbits)
