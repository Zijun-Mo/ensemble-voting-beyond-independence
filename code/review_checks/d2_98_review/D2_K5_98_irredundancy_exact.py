"""Exact rational, self-contained irredundancy witnesses for D2's 98 affine certificates.
Run in same directory as D2_K5_98_term_closed_form_candidate.py.
Proves none of these 98 terms can be deleted WITHOUT changing U_98 on
some feasible high-precision R parameter point. It does NOT prove U_98=R_det.
"""
from fractions import Fraction as F
from D2_K5_98_term_closed_form_candidate import EDGES, STATES, make_certificates

# For each of 9 S_5-orbits: (label, rational mass by binary mask, rational b by
# the fixed order EDGES=(01,02,03,04,12,13,14,23,24,34)).
WITNESSES = [
 ('dict:0', {2:'1/4',4:'1/4',8:'1/4',16:'1/4'},
  ['0']*4 + ['1/4']*6),
 ('cycle:0', {2:'1/6',4:'1/6',9:'1/6',10:'1/6',17:'1/6',20:'1/6'},
  ['0','0','1/6','1/6','0','1/6','1/6','1/6','1/6','0']),
 ('cycle:10', {5:'1/5',9:'1/5',10:'1/5',18:'1/5',20:'1/5'},
  ['0','2/5','2/5','0','0','2/5','2/5','0','2/5','0']),
 ('bp:11113', {7:'3/20',11:'3/20',13:'3/20',14:'3/20',16:'2/5'},
  ['3/10','3/10','3/10','0','3/10','3/10','0','3/10','0','0']),
 ('bp:01112', {6:'2/9',10:'2/9',12:'2/9',17:'1/3'},
  ['0','2/9','0','1/3','2/9','2/9','0','2/9','0','0']),
 ('bp:11122', {7:'1/5',8:'1/10',9:'1/10',10:'1/10',12:'1/10',16:'1/10',17:'1/10',18:'1/10',20:'1/10'},
  ['1/5','1/5','1/10','1/10','1/5','1/10','1/10','1/10','1/10','0']),
 ('bp:11111', {1:'1/20',2:'1/20',3:'3/40',4:'1/20',5:'3/40',6:'3/40',8:'1/20',9:'3/40',10:'3/40',12:'3/40',16:'1/20',17:'3/40',18:'3/40',20:'3/40',24:'3/40'},
  ['3/40']*10),
 ('rooted:0', {1:'5/26',6:'3/26',8:'3/26',10:'3/26',12:'3/26',16:'3/26',18:'3/26',20:'3/26'},
  ['0']*4+['3/26']*6),
 ('corrected:11223', {7:'1/6',11:'1/6',12:'1/4',16:'1/12',17:'1/6',18:'1/6'},
  ['1/3','1/6','1/6','1/6','1/6','1/6','1/6','1/4','0','0']),
]

CERTS = make_certificates()
assert len(CERTS) == 98
# Interior reference distribution: five independent Bernoulli(1/5) errors,
# with pair upper bounds exactly 1/25.
ref_pi = [F(1,5)**sum(int(t) for t in x)*F(4,5)**(5-sum(int(t) for t in x)) for x in STATES]
ref_b = [F(1,25)]*10
alpha=F(1,100)

def evaluate(c,q,b):
    return c.a+sum(v*t for v,t in zip(c.beta,q))+sum(v*t for v,t in zip(c.gamma,b))

def verify(label, mass, bounds):
    pi=[F(mass.get(mask, '0')) for mask in range(32)]
    b=[F(u) for u in bounds]
    assert sum(pi)==1 and all(v>=0 for v in pi)
    q0=[sum(pi[mask]*int(STATES[mask,i]) for mask in range(32)) for i in range(5)]
    assert all(b[e]>=sum(pi[mask]*int(STATES[mask,i])*int(STATES[mask,j]) for mask in range(32)) for e,(i,j) in enumerate(EDGES))
    assert all(b[e]<=min(q0[i],q0[j]) for e,(i,j) in enumerate(EDGES))
    pi=[(1-alpha)*a+alpha*t for a,t in zip(pi,ref_pi)]
    b=[(1-alpha)*a+alpha*t for a,t in zip(b,ref_b)]
    q=[sum(pi[mask]*int(STATES[mask,i]) for mask in range(32)) for i in range(5)]
    assert sum(pi)==1 and all(F(0)<t<F(1,2) for t in q)
    assert sum(q)-3*min(q)<=1
    for e,(i,j) in enumerate(EDGES):
        joint=sum(pi[mask]*int(STATES[mask,i])*int(STATES[mask,j]) for mask in range(32))
        assert joint <= b[e] <= min(q[i],q[j])
    k=next(i for i,c in enumerate(CERTS) if c.label==label)
    scores=[evaluate(c,q,b) for c in CERTS]
    gap=min(v-scores[k] for j,v in enumerate(scores) if j!=k)
    assert gap>0
    return gap

if __name__=='__main__':
    for label,mass,b in WITNESSES:
        print(label, 'strict margin =', verify(label,mass,b))
    print('Nine exact rational witnesses + relabeling cover all 98 terms.')
    print('This proves irredundancy within U98, not minimax completeness.')
