"""Compare the 19 fractional orbits listed in the K_min=9 report with the orbits found by the review's exact adjacency
decomposition (run/adj_state_K8.json), using a full S_8 canonical form (all 40320 permutations)."""
import sys, json, itertools, math, numpy as np, ast
from fractions import Fraction as F
src=open(sys.argv[1]).read()
start=src.index('ORBITS_F01_F19 = ['); end=src.index('\n]\n',start)+2
ORB=ast.literal_eval(src[start+len('ORBITS_F01_F19 = '):end])
n=8;E=list(itertools.combinations(range(n),2));idx={e:k for k,e in enumerate(E)}
perms=np.array(list(itertools.permutations(range(n))))
P=np.array([[idx[tuple(sorted((p[i],p[j])))] for (i,j) in E] for p in perms]);Pinv=np.argsort(P,axis=1)
def canon(v):
    L=1
    for x in v: L=L*x.denominator//math.gcd(L,x.denominator)
    iv=np.array([int(x*L) for x in v]);W=iv[Pinv];w=W[np.lexsort(W.T[::-1])[0]]
    return '|'.join(str(F(int(x),L)) for x in w)
S=json.load(open(sys.argv[2]))['orbits']
mine={k for k in S if any(F(x).denominator>1 for x in k.split('|'))}
theirs={canon([F(x,d) for x in t]):name for name,d,t in ORB}
print('report fractional orbits:',len(theirs),' review fractional orbits:',len(mine))
print('report orbits not found by review:',[theirs[c] for c in theirs if c not in mine])
print('review orbits not in report:',len([c for c in mine if c not in theirs]))
print('review: expanded =',sum(S[k]['expanded'] for k in S),'of',len(S),'; unexpanded:',[ (sorted(set(k.split('|'))),k.split('|').count('0')) for k in S if not S[k]['expanded']])
