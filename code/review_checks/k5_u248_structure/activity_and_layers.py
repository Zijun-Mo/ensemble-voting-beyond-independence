"""(1) Which of the 12 U248 orbits attains the minimum, over random feasible (q,b) in three regimes
(R0; R minus R0; outside R) and under two sampling designs (b at the true moments / b with random slack).
(2) On each saturation layer D_S (b_e = p_e(pi) for e in S), which orbit members are irredundant (per orbit counts).
Usage: python3 -I activity_and_layers.py <S26_bundle dir>"""
import sys, itertools, time, numpy as np
from collections import Counter
from scipy.optimize import linprog
sys.path.insert(0, sys.argv[1]); import D2_K5_U248_completeness_proof_and_verifier as U
E=list(itertools.combinations(range(5),2)); Z=np.array([[(m>>i)&1 for i in range(5)] for m in range(32)])
PAIR=np.array([[Z[m,i]*Z[m,j] for (i,j) in E] for m in range(32)])
certs=U.make_u248()
def fam(c):
    f=c.label.split(':')[0]
    if f=='bp': return 'bp'+''.join(map(str,sorted(c.weights,reverse=True)))
    if f=='cycle': return 'triangle' if sum(1 for g in c.gamma if g)==3 else 'pentagon'
    return f
FAM=[fam(c) for c in certs]; C=np.array([[float(x) for x in c.coefficients()] for c in certs])
ORDER=['dict','triangle','pentagon','bp11111','bp21110','bp22111','bp31111','rooted','corrected','anchored2','twoport','coupled_pair']
rng=np.random.default_rng(5); T=time.time()
stats={}
for design in ('moments','slack'):
    for reg in ('R0','R\\R0','outside R'): stats[(design,reg)]=Counter()
N=0
while time.time()-T<120:
    kind=rng.integers(3)
    if kind==0:
        qq=rng.uniform(0.02,0.48,5); pi=np.prod(np.where(Z==1,qq,1-qq),axis=1)
    else:
        k=rng.integers(3,10); S=rng.choice(32,k,replace=False); pi=np.zeros(32); pi[S]=rng.dirichlet(np.ones(k))
        if kind==2: qq=rng.uniform(0.05,0.45,5); pi=0.5*pi+0.5*np.prod(np.where(Z==1,qq,1-qq),axis=1)
    q=Z.T@pi
    if np.any(q>=0.5): continue
    reg='R0' if q.sum()<=1 else ('R\\R0' if q.sum()-3*q.min()<=1 else 'outside R')
    p=PAIR.T@pi; mq=np.array([min(q[i],q[j]) for i,j in E])
    for design in ('moments','slack'):
        b=p if design=='moments' else p+rng.random(10)*(rng.random(10)<0.5)*(mq-p)
        v=C[:,0]+C[:,1:6]@q+C[:,6:]@b; m=v.min(); win={FAM[j] for j in np.where(v<=m+1e-12)[0]}
        stats[(design,reg)][min(win,key=ORDER.index) if len(win)>1 else next(iter(win))]+=1
    N+=1
print(f'{N} random feasible parameter points (independent / sparse / mixed joint laws)')
for (design,reg),c in stats.items():
    tot=sum(c.values())
    print(f'  b={design:7s} {reg:9s} n={tot:5d}: '+', '.join(f'{f} {100*c[f]/tot:.1f}%' for f in ORDER if c[f]))
# (2) layers
idx={e:k for k,e in enumerate(E)}
layers={'all (no saturation)':set(),'disjoint M={01,23}':set(range(10))-{idx[(0,1)],idx[(2,3)]},'adjacent M={01,02}':set(range(10))-{idx[(0,1)],idx[(0,2)]},
        'single M={01}':set(range(10))-{idx[(0,1)]},'all saturated':set(range(10))}
LP=C[:,[0]]+C[:,1:6]@Z.T; LB=C[:,6:]
for name,S in layers.items():
    A=[];bb=[];Aeq=[np.r_[np.ones(32),np.zeros(11)]];beq=[1.0]
    for e,(i,j) in enumerate(E):
        z=np.zeros(10); z[e]=-1
        (Aeq if e in S else A).append(np.r_[PAIR[:,e],z,0]); (beq if e in S else bb).append(0.0)
        for k in (i,j):
            z=np.zeros(10); z[e]=1; A.append(np.r_[-Z[:,k],z,0]); bb.append(0.0)
    for k in range(5): A.append(np.r_[Z[:,k],np.zeros(10),0]); bb.append(0.5)
    A=np.array(A); bb=np.array(bb); keep=Counter()
    for j in range(248):
        o=[i for i in range(248) if i!=j]
        Rr=np.c_[-(LP[o]-LP[j]),-(LB[o]-LB[j]),np.ones(len(o))]
        r=linprog(np.r_[np.zeros(42),-1.0],A_ub=np.vstack([A,Rr]),b_ub=np.r_[bb,np.zeros(len(o))],A_eq=np.array(Aeq),b_eq=beq,bounds=[(0,None)]*42+[(None,None)],method='highs')
        if -r.fun>1e-7: keep[FAM[j]]+=1
    tot=Counter(FAM)
    print(f'  layer {name:22s}: {sum(keep.values()):3d} irredundant; dropped orbits/members: '+', '.join(f'{f} {tot[f]-keep[f]}/{tot[f]}' for f in ORDER if keep[f]<tot[f]))
