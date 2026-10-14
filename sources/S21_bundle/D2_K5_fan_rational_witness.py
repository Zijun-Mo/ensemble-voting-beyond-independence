"""Exact rational, constructive minimax witness for the five-vertex fan.

Checks an instance outside the earlier triangle-surplus condition and
constructs a common adversarial distribution (nu + rho) on all 32 states.
Only Python's standard-library fractions; no numerical optimizer required.
"""
from fractions import Fraction as F
from itertools import combinations

E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(2,3)]
TRI=[(0,1,2),(0,1,4),(0,2,3)]
q=[F(3,10),F(8,25),F(23,50),F(23,50),F(3,10)]
eps=F(1,5)

# A distribution witnessing feasibility of the original marginal/cap data.
source=[]
for m in range(32):
    only={0:1-q[0]-q[1]-q[4],1:q[0],2:q[1],16:q[4]}.get(m & 19,F(0))
    for i in (2,3): only*=q[i] if (m>>i)&1 else 1-q[i]
    independent=F(1)
    for i in range(5): independent*=q[i] if (m>>i)&1 else 1-q[i]
    source.append((1-eps)*only+eps*independent)
b=[sum(source[m] for m in range(32) if (m>>i)&1 and (m>>j)&1)
   for i,j in E]
assert sum(source)==1
assert all(sum(source[m] for m in range(32) if (m>>i)&1)==q[i]
           for i in range(5))

# Rational t,u from the projected 8-variable LP (verified below by
# explicitly constructing its matching joint law, not trusting the LP).
t=F(141,2500)
u=[F(12,625),F(93,2500),F(0),F(9,500),
   F(0),F(12,625),F(12,625)]
h=1-2*t
r=[qi-t for qi in q]
d=[be-ue for be,ue in zip(b,u)]

def marg(law,ids):
    mask=sum(1<<i for i in ids)
    return {z:sum(mass for m,mass in law.items() if (m&mask)==z)
            for z in range(32) if (z&mask)==z}

def glue(law_a,nodes_a,law_b,nodes_b):
    """Clique-tree gluing of two unnormalized laws sharing a separator."""
    overlap=set(nodes_a)&set(nodes_b)
    common=sum(1<<i for i in overlap)
    sep=marg(law_a,overlap)
    assert sep==marg(law_b,overlap)
    joined={}
    for ma,pa in law_a.items():
        for mb,pb in law_b.items():
            if (ma&common)!=(mb&common): continue
            den=sep[ma&common]
            if den==0: assert pa==0 or pb==0;continue
            key=ma|mb
            joined[key]=joined.get(key,F(0))+pa*pb/den
    return joined

def symmetric_triangle(nodes):
    ids=[E.index(edge) for edge in combinations(nodes,2)]
    uu=[u[k] for k in ids]
    total=sum(uu)
    both=(total-t)/2
    out={}
    def push(mask,prob):out[mask]=out.get(mask,F(0))+prob
    full=sum(1<<i for i in nodes)
    push(0,both);push(full,both)
    for k,(i,j) in zip(ids,combinations(nodes,2)):
        mass=(2*u[k]-total+t)/2
        on=(1<<i)|(1<<j)
        push(on,mass);push(full^on,mass)
    assert all(z>=0 for z in out.values())
    assert sum(out.values())==2*t
    return out

def central_residual():
    v01=v02=F(0)
    out={}
    for val in (0,1):
        universe=r[0] if val else h-r[0]
        s1=v01 if val else r[1]-v01
        s2=v02 if val else r[2]-v02
        both=max(F(0),s1+s2-universe)
        mass={(0,0):universe-s1-s2+both,
              (1,0):s1-both,(0,1):s2-both,(1,1):both}
        for (v1,v2),p in mass.items():
            out[val|(v1<<1)|(v2<<2)]=p
    assert all(z>=0 for z in out.values())
    return out

def leaf_residual(i,j,k,shared_v,cap_ik,cap_jk):
    # Greedy: put X_k=1 first into the 00 cell (free), then into cells
    # 10 / 01 (one budget each), and finally 11 (both budgets).
    p={(0,0):h-r[i]-r[j]+shared_v,(1,0):r[i]-shared_v,
       (0,1):r[j]-shared_v,(1,1):shared_v}
    assert all(v>=0 for v in p.values())
    alloc={(0,0):p[0,0]}
    alloc[1,0]=min(p[1,0],cap_ik)
    alloc[0,1]=min(p[0,1],cap_jk)
    alloc[1,1]=min(p[1,1],cap_ik-alloc[1,0],cap_jk-alloc[0,1])
    assert min(alloc.values())>=0
    max_mass=sum(alloc.values())
    assert max_mass>=r[k]
    fraction=r[k]/max_mass if max_mass else F(0)
    out={}
    for (xi,xj),size in p.items():
        yes=alloc[xi,xj]*fraction
        base=(xi<<i)|(xj<<j)
        out[base]=size-yes
        out[base|(1<<k)]=yes
    assert all(z>=0 for z in out.values())
    return out

nu=glue(symmetric_triangle(TRI[0]),TRI[0],
        symmetric_triangle(TRI[1]),TRI[1])
nu=glue(nu,[0,1,2,4],symmetric_triangle(TRI[2]),TRI[2])
rho=glue(central_residual(),TRI[0],
         leaf_residual(0,1,4,F(0),d[3],d[5]),TRI[1])
rho=glue(rho,[0,1,2,4],
         leaf_residual(0,2,3,F(0),d[2],d[6]),TRI[2])
pi=[nu.get(m,F(0))+rho.get(m,F(0)) for m in range(32)]
assert all(x>=0 for x in pi)
assert sum(pi)==1 and sum(nu.values())==2*t
assert all(nu.get(m,F(0))==nu.get(31-m,F(0)) for m in range(32))
assert all(sum(pi[m] for m in range(32) if m>>i&1)==q[i] for i in range(5))
assert all(sum(pi[m] for m in range(32) if m>>i&1 and m>>j&1)<=b[k]
           for k,(i,j) in enumerate(E))
cycles=[[0,1,4],[0,3,5],[1,2,6],[2,3,4,5,6]]
candidate=min([*q,*[sum(b[k] for k in cycle) for cycle in cycles]])
surplus=max(sum(q[i] for i in tri)-1 for tri in TRI)
assert candidate==t and t<surplus
assert sum(q)-3*min(q)<1 and sum(q)>1
print('RATIONAL THEOREM WITNESS')
print('q:',[str(z) for z in q])
print('b:',[str(z) for z in b])
print('candidate and achieved minimax:',t)
print('old sufficient-condition RHS:',surplus)
print('nu mass:',sum(nu.values()),'; full distribution support:',sum(z>0 for z in pi))
print('ALL 32 STATE PROBABILITIES, MARGINALS AND EDGE CAPS CHECK EXACTLY')
