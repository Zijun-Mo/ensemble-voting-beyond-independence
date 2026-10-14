"""Structure of U248 (result 8): for each of the 12 S5-orbits give a representative certificate, its intercept a,
the zero-intercept part Q0 = Q + a*(N-1) (N = number of errors), its interaction graph, and ALL self-dual threshold rules
(up to S5, by weight type, with signs) that the certificate dominates (Q >= F on all 32 states).
Usage: python3 -I orbit_table.py <S26_bundle dir>"""
import sys, itertools, json
import networkx as nx
from fractions import Fraction as F
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U
E=list(itertools.combinations(range(5),2)); Z=[[(m>>i)&1 for i in range(5)] for m in range(32)]
CAN=[(1,0,0,0,0),(1,1,1,0,0),(1,1,1,1,1),(2,1,1,1,0),(2,2,1,1,1),(3,1,1,1,1),(3,2,2,1,1)]
rules={}
for w0 in CAN:
    for w in set(itertools.permutations(w0)):
        for sg in itertools.product((1,-1),repeat=5):
            ww=tuple(a*b for a,b in zip(w,sg)); f=tuple(int(sum(ww[i]*(2*z[i]-1) for i in range(5))>0) for z in Z)
            rules.setdefault(f,ww)
assert len(rules)==1882
def Qval(c,z): return c[0]+sum(c[1+i]*z[i] for i in range(5))+sum(c[6+e]*z[i]*z[j] for e,(i,j) in enumerate(E))
def gname(G):
    n,m=G.number_of_nodes(),G.number_of_edges()
    if m==0: return 'empty'
    cands={'K3':nx.complete_graph(3),'C5':nx.cycle_graph(5),'K4':nx.complete_graph(4),'K5':nx.complete_graph(5),
           'K5-e':nx.complete_graph(5),'gem (fan)':nx.Graph([(0,1),(1,2),(2,3),(4,0),(4,1),(4,2),(4,3)])}
    cands['K5-e'].remove_edge(0,1)
    for k,H in cands.items():
        if nx.is_isomorphic(G,H): return k
    return f'{n} vertices, {m} edges'
seen={}; out=[]
P=list(itertools.permutations(range(5)))
for c in U.make_u248():
    fam=c.label.split(':')[0]
    if fam=='bp': fam='bp'+''.join(map(str,sorted(c.weights,reverse=True)))
    if fam=='cycle': fam='triangle' if sum(1 for g in c.gamma if g)==3 else 'pentagon'
    if fam in seen: seen[fam]+=1; continue
    seen[fam]=1; v=c.coefficients()
    a=v[0]; beta=v[1:6]; gam=v[6:]
    beta0=[b+a for b in beta]
    G=nx.Graph([E[e] for e in range(10) if gam[e]!=0])
    dom=sorted({tuple(sorted(map(abs,w),reverse=True))+('neg' if min(w)<0 else '',) for f,w in rules.items() if all(Qval(v,Z[m])>=f[m] for m in range(32))})
    tight=[m for m in range(32) if Qval(v,Z[m])==int(2*sum(c.weights[i]*Z[m][i] for i in range(5))>sum(c.weights))]
    out.append(dict(family=fam,label=c.label,a=str(a),beta=[str(x) for x in beta],gamma={f'{i}{j}':str(gam[e]) for e,(i,j) in enumerate(E) if gam[e]!=0},
                    Q0_beta=[str(x) for x in beta0],graph=gname(G),rule_weights=list(c.weights),dominated_rule_types=[list(d) for d in dom],
                    l1=str(sum(abs(x) for x in beta)+sum(gam))))
for fam in out: fam['orbit_size']=seen[fam['family']]
for o in out:
    print(f"{o['family']:13s} size {o['orbit_size']:2d}  a={o['a']:4s} beta={o['beta']}  Q0_beta={o['Q0_beta']}  graph={o['graph']}\n"
          f"              gamma={o['gamma']}\n              rules dominated: {o['dominated_rule_types']}  l1={o['l1']}")
json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else '/dev/null','w'),indent=1)
