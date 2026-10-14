"""D2, K=5: an explicit 98-affine-term robust voting certificate family.

THIS IS A CANDIDATE FOR AN EXACT MINIMAX CLOSED FORM, NOT A PROVED EQUALITY.
Each of its 98 terms is rigorously a valid (possibly non-sharp) worst-case
error upper bound, under feasible marginal / pairwise joint-error caps.

The terms are: 5 dictators + 22 odd cycles + 36 weighted-BP polynomials
+ 5 rooted fractional polynomials + 30 corrected heavy-triplet polynomials.

No LP is required to evaluate this explicit bound. Validity checks evaluate
every polynomial on all 32 binary error patterns.

Requires only Python and numpy (for checks). Compatible with Python>=3.10.
"""
from __future__ import annotations
import itertools
from dataclasses import dataclass
from fractions import Fraction
import numpy as np

K=5
EDGES=tuple(itertools.combinations(range(K),2))
STATES=np.array([[(mask>>i)&1 for i in range(K)] for mask in range(32)],dtype=int)
PAIR=np.array([STATES[:,i]*STATES[:,j] for i,j in EDGES]).T
MOMENT=np.column_stack([np.ones(32,dtype=int),STATES,PAIR])

@dataclass(frozen=True)
class Certificate:
    label: str
    a: Fraction
    beta: tuple[Fraction,...]
    gamma: tuple[Fraction,...]
    weights: tuple[int,...]

    def coefficients(self): return (self.a,)+self.beta+self.gamma
    def bound(self,q,b):
        return self.a + sum(v*Fraction(str(x)) for v,x in zip(self.beta,q)) + \
               sum(v*Fraction(str(x)) for v,x in zip(self.gamma,b))
    def check(self):
        assert len(self.beta)==5 and len(self.gamma)==10
        assert all(z>=0 for z in self.gamma)
        assert sum(self.weights)%2==1 and sum(self.weights)>0
        for x in STATES:
            L=self.a + sum(c*int(t) for c,t in zip(self.beta,x)) + \
                       sum(c*int(x[i])*int(x[j]) for c,(i,j) in zip(self.gamma,EDGES))
            F=int(2*sum(int(v)*int(t) for v,t in zip(self.weights,x))>sum(self.weights))
            assert L>=F,(self.label,tuple(x),L,F)


def odd_cycles():
    cc=[]
    for vertices in itertools.combinations(range(K),3):
        cc.append(tuple(EDGES.index(e) for e in itertools.combinations(vertices,2)))
    for rest in itertools.permutations(range(1,K)):
        cyc=(0,)+rest
        if cyc[1]>cyc[-1]:continue
        cc.append(tuple(EDGES.index(tuple(sorted((cyc[i],cyc[(i+1)%K])))) for i in range(K)))
    assert len(cc)==22
    return cc


def weighted_bp(weights,label):
    w=tuple(weights);total=sum(w)
    assert total>=3 and total%2==1
    r=(total+1)//2
    a=Fraction((r-2)*(r-1),2)
    beta=tuple(Fraction(ww*(ww-2*r+3),2) for ww in w)
    gamma=tuple(Fraction(w[i]*w[j]) for i,j in EDGES)
    return Certificate(label,a,beta,gamma,w)


def rooted(i):
    w=tuple(3 if j==i else 1 for j in range(K))
    a=Fraction(1,3)
    beta=(Fraction(-1,3),)*K
    gamma=tuple(Fraction(4 if i in (j,k) else 1,3) for j,k in EDGES)
    return Certificate('rooted:%d'%i,a,beta,gamma,w)


def corrected_core(weights):
    w=tuple(weights)
    assert sorted(w)==[1,1,2,2,3]
    root=weighted_bp(w,'unused')
    H={i for i,wi in enumerate(w) if wi>=2}
    beta=tuple(root.beta[i]+(1 if i in H else 0) for i in range(K))
    gamma=tuple(root.gamma[e]-(1 if i in H and j in H else 0)
                for e,(i,j) in enumerate(EDGES))
    return Certificate('corrected:%s'%(''.join(map(str,w))),root.a-1,beta,gamma,w)


def make_certificates():
    out=[]
    for i in range(K):
        w=tuple(int(j==i) for j in range(K))
        out.append(Certificate('dict:%d'%i,Fraction(0),tuple(Fraction(v) for v in w),(Fraction(0),)*10,w))
    for ix,C in enumerate(odd_cycles()):
        # Odd cycle votes are a committee of 3 or 5, irrespective of which cyclic edges bound its error.
        vertices={i for e in C for i in EDGES[e]}
        w=tuple(int(i in vertices) for i in range(K))
        gamma=tuple(Fraction(int(e in C)) for e in range(10))
        out.append(Certificate('cycle:%d'%ix,Fraction(0),(Fraction(0),)*5,gamma,w))
    for canonical in ((3,1,1,1,1),(2,1,1,1,0),(2,2,1,1,1),(1,1,1,1,1)):
        for w in sorted(set(itertools.permutations(canonical))):
            out.append(weighted_bp(w,'bp:%s'%(''.join(map(str,w)))))
    for i in range(K):out.append(rooted(i))
    for w in sorted(set(itertools.permutations((3,2,2,1,1)))):
        out.append(corrected_core(w))
    assert len(out)==98,len(out)
    assert len(set(c.coefficients() for c in out))==98
    return out


def bound(q,b):
    candidates=[(c.bound(q,b),c.label) for c in make_certificates()]
    return min(candidates)

if __name__=='__main__':
    cert=make_certificates()
    for c in cert:c.check()
    print('Verified all 98 polynomial majorants exactly on all 32 patterns.')
    q=[Fraction(9,26)]*4+[Fraction(5,26)]
    b=[Fraction(3,26) if 4 not in e else Fraction(0) for e in EDGES]
    v,label=bound(q,b)
    assert v==Fraction(1,26)
    print('Targeted example: bound=',v,'label=',label)
