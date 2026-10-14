"""Independent exact verifier (written by the reviewer) for rational certificates of B_graph(u) >= 1 on K_8.

Certificate = binary decision tree. Internal node {"cut": S} (S = bitmask of a vertex subset of {0..7}) has children
"A" (assumes y(E(S)) >= 1) and "D" (assumes y(E(V\\S)) >= 1). Leaf {"leaf": [[j, "lambda_j"], ...]} where j indexes the
constraints on the root-to-leaf path (0 = root's). Leaf checks (pure Fractions):
    lambda_j >= 0,  sum_j lambda_j >= 1,  and for every edge e:  sum_j lambda_j [e in side_j] <= u_e.
Soundness: for y >= 0 with g(y) >= 1, every node's cut has max(y(E(S)), y(E(V\\S))) >= 1, so y follows some root-to-leaf
path; at that leaf u.y >= sum_j lambda_j y(side_j) >= sum_j lambda_j >= 1. Hence B_graph(u) >= 1; together with
B_graph <= c* = 1 this gives B_graph(u) = 1.
Usage: python3 -I verify_B_certs.py <orbit file (S27 script)> <cert_dir>   (checks every <cert_dir>/*.json)"""
import sys, json, os, ast, itertools
from fractions import Fraction as F
E=list(itertools.combinations(range(8),2))
src=open(sys.argv[1]).read(); a=src.index('ORBITS_F01_F19 = ['); b=src.index('\n]\n',a)+2
ORB={name:[F(x,d) for x in t] for name,d,t in ast.literal_eval(src[a+len('ORBITS_F01_F19 = '):b])}
def side(S,which):   # edge indicator of E(S) (which='A') or E(V\S) (which='D')
    inS=[(S>>v)&1 for v in range(8)]
    want=1 if which=='A' else 0
    return [int(inS[i]==want and inS[j]==want) for i,j in E]
def check(node,path,u,stats):
    if 'leaf' in node:
        lam=[(int(j),F(l)) for j,l in node['leaf']]
        assert all(0<=j<len(path) for j,_ in lam), 'bad index'
        assert all(l>=0 for _,l in lam), 'negative multiplier'
        assert sum(l for _,l in lam)>=1, 'sum of multipliers < 1'
        for e in range(28):
            assert sum(l*path[j][e] for j,l in lam)<=u[e], f'edge {E[e]} over budget'
        stats['leaves']+=1; stats['depth']=max(stats['depth'],len(path))
        stats['maxden']=max([stats['maxden']]+[l.denominator for _,l in lam])
        return
    S=int(node['cut']); assert 0<S<255
    check(node['A'],path+[side(S,'A')],u,stats)
    check(node['D'],path+[side(S,'D')],u,stats)
ok=0
for fn in sorted(os.listdir(sys.argv[2])):
    if not fn.endswith('.json'): continue
    c=json.load(open(os.path.join(sys.argv[2],fn)))
    name=c['orbit']; u=[F(x) for x in c['u']]
    assert u==ORB[name], f'{name}: certificate vector differs from the listed orbit representative'
    st={'leaves':0,'depth':0,'maxden':1}
    check(c['tree'],[],u,st); ok+=1
    print(f"{name}: B_graph(u) >= 1 certified exactly; leaves {st['leaves']}, max depth {st['depth']}, max denominator {st['maxden']}")
print(f'{ok} certificates verified; orbits without certificate: {sorted(set(ORB)-{json.load(open(os.path.join(sys.argv[2],f)))["orbit"] for f in os.listdir(sys.argv[2]) if f.endswith(".json")})}')
