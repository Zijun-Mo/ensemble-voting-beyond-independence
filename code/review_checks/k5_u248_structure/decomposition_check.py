"""Intercept decomposition of U248: every certificate Q = Q0 - a*(N-1), where N = number of errors and Q0 has zero
intercept. Checks EXACTLY that each Q0 is a legal zero-intercept certificate for the same rule (so L_Q0 >= U0 =
min{q_min, c*} by result 3), and that its linear coefficients are nonnegative. Corollary checked numerically:
U0 - 5*(sum q - 1)_+ <= R*_det = U248 <= U0. Also prints the 12 exposure witnesses of [S26] (q, sum q, winning margin).
Usage: python3 -I decomposition_check.py <S26_bundle dir>"""
import sys, itertools, numpy as np
from fractions import Fraction as F
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U
E=list(itertools.combinations(range(5),2)); Z=[[(m>>i)&1 for i in range(5)] for m in range(32)]
certs=U.make_u248(); amax=0
for c in certs:
    v=c.coefficients(); a=v[0]
    q0=[F(0)]+[v[1+i]+a for i in range(5)]+list(v[6:])
    w=c.weights
    for z in Z:
        Fx=int(2*sum(w[i]*z[i] for i in range(5))>sum(w))
        Q0=sum(q0[1+i]*z[i] for i in range(5))+sum(q0[6+e]*z[i]*z[j] for e,(i,j) in enumerate(E))
        assert Q0>=Fx, (c.label,z)
    assert all(x>=0 for x in q0[1:]) and a>=0
    amax=max(amax,a)
print('all 248: Q0 = Q + a(N-1) is a legal zero-intercept certificate for the same rule, coefficients >= 0; max a =',amax)
# corollary
C=np.array([[float(x) for x in c.coefficients()] for c in certs]); zero=[j for j,c in enumerate(certs) if c.a==0]
rng=np.random.default_rng(9); Zn=np.array(Z); worst=0; n=0
while n<200000:
    k=rng.integers(3,10); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k))
    q=Zn.T@pi
    if np.any(q>=0.5): continue
    p=np.array([pi@(Zn[:,i]*Zn[:,j]) for i,j in E]); b=p+rng.random(10)*np.array([min(q[i],q[j]) for i,j in E]-p)*(rng.random(10)<0.4)
    v=C[:,0]+C[:,1:6]@q+C[:,6:]@b; U0=v[zero].min(); R=v.min(); tau=max(q.sum()-1,0)
    assert R<=U0+1e-12 and R>=U0-5*tau-1e-12
    if tau>0: worst=max(worst,(U0-R)/tau)
    n+=1
print(f'{n} random points: U0 - 5 tau <= R* <= U0 always; largest observed (U0-R*)/tau = {worst:.3f}')
delta=F(1,1000); anc=[F(1,5)**sum(z)*F(4,5)**(5-sum(z)) for z in Z]
for idx,label,sz,pnz,bnz in U.SPARSE_12_ORBIT_WITNESSES:
    pm=[(1-delta)*F(pnz.get(m,'0'))+delta*anc[m] for m in range(32)]
    q=[sum(pm[m]*Z[m][i] for m in range(32)) for i in range(5)]
    vals=sorted(float(c.bound(q,[(1-delta)*F(bnz.get(e,'0'))+delta*(F(1,25)+F(1,100)) for e in range(10)])) for c in certs)
    print(f'  witness {label:40s} q={[round(float(x),3) for x in q]} sum q={float(sum(q)):.3f}  value {vals[0]:.4f}, next {vals[1]:.4f}')
