"""Exact construction for a C5 with one chord (01,12,23,34,04,02).

Theorem: if 0<=q_i<1/2, non-chord edges are uninformative,
 t=min(min(q),b01+b12+b02,b01+b12+b23+b34+b04),
and t>=q0+q1+q2-1, then R_det^*=t.
Constructs the common adversarial distribution using Fraction only.
"""
from fractions import Fraction as F
from itertools import combinations
from collections import defaultdict

V=range(5)
E=list(combinations(V,2))
H={(0,1),(1,2),(0,2),(2,3),(3,4),(0,4)}
C3={(0,1),(1,2),(0,2)}
C5={(0,1),(1,2),(2,3),(3,4),(0,4)}

def mask_for_defects(defects):
    target=frozenset(defects)
    fits=[]
    for mask in range(32):
        agree=frozenset(e for e in H if ((mask>>e[0])&1)==((mask>>e[1])&1))
        if agree==target: fits.append(mask)
    assert len(fits)==2 and fits[0]^fits[1]==31,(target,fits)
    return fits

def construct(q,b):
    """q sequence, b dictionary indexed by canonical edges, exact rationals."""
    q=[F(v) for v in q]
    b={e:F(v) for e,v in b.items()}
    assert 0<=min(q) and max(q)<F(1,2)
    assert all(0<=b[e]<=min(q[e[0]],q[e[1]]) for e in E)
    assert all(b[e]==min(q[e[0]],q[e[1]]) for e in E if e not in H)
    t=min(min(q),sum(b[e] for e in C3),sum(b[e] for e in C5))
    assert t>=q[0]+q[1]+q[2]-1,'Additional triangle condition is required'

    # Symmetric pairs: each pair carries 2*lambda total mass, and
    # contributes lambda to E[X_i X_j] exactly when edge ij is a defect.
    r=t
    first=min(r,b[(0,1)]);r-=first
    second=min(r,b[(1,2)]);r-=second
    assert r<=b[(0,2)]
    third=min(r,b[(2,3)]);r-=third
    fourth=min(r,b[(3,4)]);r-=fourth
    fifth=r
    assert 0<=fifth<=b[(0,4)]
    slots=[({(0,1)},first),({(1,2)},second),
           ({(0,2),(2,3)},third),({(0,2),(3,4)},fourth),
           ({(0,2),(0,4)},fifth)]
    assert sum(mass for _,mass in slots)==t
    pi=defaultdict(F)
    for defects,mass in slots:
        if mass:
            masks=mask_for_defects(defects)
            for m in masks:pi[m]+=mass

    # Residual subdistribution: every supported set is H-independent.
    # Separator {0,2}: the triangle {0,1,2} and square {0,2,3,4}
    # meet in the edge 02.  Work with UNNORMALIZED residual marginals.
    z=1-2*t
    d=[v-t for v in q]
    assert all(v>=0 for v in d)
    assert d[0]+d[1]+d[2]<=z
    assert all(d[i]+d[j]<=z for i,j in H)

    # First allocate a measure on (X0,X2,X3,X4), no adjacent ones.
    base=defaultdict(F)
    def add(mask,mass):
        assert mass>=0
        if mass:base[mask]+=mass
    a=min(d[3],d[0]); add((1<<0)|(1<<3),a);add(1<<0,d[0]-a)
    c=min(d[4],d[2]); add((1<<2)|(1<<4),c);add(1<<2,d[2]-c)
    remain=z-d[0]-d[2]
    extra3=d[3]-a;extra4=d[4]-c
    assert remain>=extra3+extra4>=0
    add(1<<3,extra3);add(1<<4,extra4)
    add(0,remain-extra3-extra4)
    assert sum(base.values())==z

    # Add x1 only where x0=x2=0, preserving all graph independence.
    none_mass=sum(p for mask,p in base.items() if not(mask&((1<<0)|(1<<2))))
    assert none_mass==remain and d[1]<=remain
    for mask,p in base.items():
        if mask&((1<<0)|(1<<2)):
            pi[mask]+=p
        elif remain:
            mass1=p*d[1]/remain
            pi[mask|(1<<1)]+=mass1
            pi[mask]+=p-mass1
        else:
            assert d[1]==0
            pi[mask]+=p

    assert sum(pi.values())==1
    q_actual=[sum(v*((mask>>i)&1) for mask,v in pi.items()) for i in V]
    assert q_actual==q,(q_actual,q)
    pair_actual={e:sum(v*((mask>>e[0])&1)*((mask>>e[1])&1)
                       for mask,v in pi.items()) for e in E}
    assert all(pair_actual[e]<=b[e] for e in E),(pair_actual,b)
    assert sum(pi.values())==1
    # The t symmetric pair mass simultaneously yields failure >=t
    # for EVERY deterministic self-dual rule.
    return t,dict(sorted(pi.items())),pair_actual

def examples():
    # A genuinely heterogeneous R\R0 example with t < min(q_i)
    q=[F(28,100),F(25,100),F(30,100),F(28,100),F(25,100)]
    caps={(0,1):F(2,100),(1,2):F(3,100),(0,2):F(8,100),
          (2,3):F(5,100),(3,4):F(4,100),(0,4):F(2,100)}
    b={e:caps[e] if e in H else min(q[e[0]],q[e[1]]) for e in E}
    t,pi,actual=construct(q,b)
    assert t==F(13,100) and sum(q)>1
    assert sum(q)-3*min(q)<=1
    assert len(pi)<=32
    print('EXAMPLE 1: q=',[str(v) for v in q],'t=',t,
          'number of nonzero atoms=',len(pi))
    print('  distribution:',{m:str(v) for m,v in pi.items() if v})

    # Triangle-overload case, where q0+q1+q2>1 but theorem still works.
    q=[F(42,100),F(40,100),F(40,100),F(40,100),F(41,100)]
    caps={(0,1):F(8,100),(1,2):F(8,100),(0,2):F(10,100),
          (2,3):F(4,100),(3,4):F(4,100),(0,4):F(5,100)}
    b={e:caps[e] if e in H else min(q[e[0]],q[e[1]]) for e in E}
    t,pi,actual=construct(q,b)
    assert t==F(26,100) and t>=q[0]+q[1]+q[2]-1
    assert sum(q)-3*min(q)<=1
    print('EXAMPLE 2: q=',[str(v) for v in q],'t=',t,
          'triangle overload=',sum(q[i] for i in (0,1,2))-1)
    # Feasible data outside the sufficient condition: the theorem does NOT
    # claim the formula beyond that condition.  This exact witness prevents
    # conflating a sufficient region with full feasibility.
    pi_bad={5:F(3,20),9:F(3,20),10:F(3,20),18:F(3,20),20:F(3,20),
            11:F(1,20),13:F(1,20),21:F(1,20),22:F(1,20),26:F(1,20)}
    q_bad=[F(9,20)]*5
    b_bad={e:(F(1,20) if e in C5 else F(2,5) if e==(0,2)
              else F(9,20)) for e in E}
    assert sum(pi_bad.values())==1
    assert all(sum(p*((m>>i)&1) for m,p in pi_bad.items())==q_bad[i]
               for i in V)
    assert all(sum(p*((m>>i)&1)*((m>>j)&1) for m,p in pi_bad.items())<=b_bad[i,j]
               for i,j in E)
    t_bad=min(min(q_bad),sum(b_bad[e] for e in C3),sum(b_bad[e] for e in C5))
    assert t_bad==F(1,4)<sum(q_bad[i] for i in (0,1,2))-1
    print('BOUNDARY CHECK: feasible but outside sufficient condition, t=',t_bad)
    print('All construction checks PASSED.')

if __name__=='__main__':examples()
