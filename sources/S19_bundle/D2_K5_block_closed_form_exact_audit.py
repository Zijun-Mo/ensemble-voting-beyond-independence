"""Exact census and example for the proven K5 block/cycle closed-form region.
The full algebraic proof requires the K4 mixed 24-vertex exact C++ audit and
block gluing theorem, NOT empirical testing. Uses only stdlib + networkx.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from collections import Counter, defaultdict
import networkx as nx
V=range(5);E=list(combinations(V,2))
count=Counter();orbits=defaultdict(set)
for mask in range(1<<10):
    G=nx.Graph();G.add_nodes_from(V);G.add_edges_from(e for j,e in enumerate(E) if mask>>j&1)
    blocks=list(nx.biconnected_components(G))
    has_spanning_block=any(len(B)==5 for B in blocks)
    is_C5=(G.number_of_edges()==5 and all(G.degree(i)==2 for i in V) and nx.is_connected(G))
    category=('solved_C5' if is_C5 else
              'solved_blocks_at_most_4' if not has_spanning_block else
              'solved_bipartite_block5' if nx.is_bipartite(G) else
              'remaining_biconnected')
    count[category]+=1
    edges=list(G.edges)
    canonical=min(sum(1<<E.index(tuple(sorted((p[i],p[j])))) for i,j in edges)
                  for p in permutations(V))
    orbits[category].add(canonical)
assert count=={'solved_blocks_at_most_4':786,'solved_C5':12,'solved_bipartite_block5':10,'remaining_biconnected':216},count
assert {k:len(v) for k,v in orbits.items()}=={
    'solved_blocks_at_most_4':24, 'solved_C5':1,
    'solved_bipartite_block5':1, 'remaining_biconnected':8}
print('ORBIT COUNTS', {k:len(v) for k,v in orbits.items()})
print('PASS graph classes:',dict(count),'; solved',count['solved_blocks_at_most_4']+count['solved_C5']+count['solved_bipartite_block5'],'/1024')
print('PASS unresolved graph isomorphism classes = 8')
# Explicit rational example: effective H=K4 + pendant edge 04.
pi={1:F(1,5),1|16:F(1,5)}
for i,j in combinations((1,2,3),2):
    m=(1<<i)|(1<<j)
    pi[m]=F(4,25)
    pi[m|16]=F(1,25)
assert sum(pi.values())==1
q={i:sum(p*((m>>i)&1) for m,p in pi.items()) for i in V}
actual={e:sum(p*((m>>e[0]&1)*(m>>e[1]&1)) for m,p in pi.items()) for e in E}
b={e:actual[e] for e in E}
for j in (1,2,3):b[tuple(sorted((4,j)))]=min(q[4],q[j])
assert all(actual[e]<=b[e]<=min(q[e[0]],q[e[1]]) for e in E)
assert q=={0:F(2,5),1:F(2,5),2:F(2,5),3:F(2,5),4:F(8,25)}
H=nx.Graph();H.add_nodes_from(V);H.add_edges_from(e for e in E if b[e]<min(q[e[0]],q[e[1]]))
assert set(H.edges)=={(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(2,3)}
assert sum(q.values())-3*min(q.values())==F(24,25)
# Anchor (2,1,1,1) on core {0,1,2,3}: failure x0=1 and someone else error,
# or x0=0 and all three others err. Neither occurs on support of pi.
for mask in pi:
    N=2*((mask>>0)&1)+sum((mask>>i)&1 for i in (1,2,3))
    assert N<3
H_root=1-sum(q[i] for i in range(4)) + sum(b[e] for e in combinations((1,2,3),2)) + 2*sum(b[(0,i)] for i in (1,2,3))
assert H_root==0
triangles=[sum(b[e] for e in combinations(S,2)) for S in combinations(V,3) if all(H.has_edge(*e) for e in combinations(S,2))]
assert min(triangles)==F(1,5)
print('PASS rational K4+pendant example: q=',list(q.values()),'qmin=',min(q.values()),
      'cheapest triangle=',min(triangles),'anchor=',H_root,'R condition=24/25')
