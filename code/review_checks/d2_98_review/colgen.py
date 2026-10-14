"""Find (q,b) in R where some valid quadratic certificate beats the current catalog; verify exactly in rationals;
store new certificates (as S5 orbits). Usage: python3 colgen.py seed seconds"""
import numpy as np, time, sys, json, itertools, os
from fractions import Fraction as Fr
from exact_core import *
from D2_K5_98_term_closed_form_candidate import make_certificates
PERMS=list(itertools.permutations(range(5)))
def coef_perm(c,p):  # c=(a,beta5,gamma10) ints/fractions; relabel i->p[i]
    a=c[0]; beta=[0]*5; gam=[0]*10
    for i in range(5): beta[p[i]]=c[1+i]
    for e,(i,j) in enumerate(EDGES): gam[EDGES.index(tuple(sorted((p[i],p[j]))))]=c[6+e]
    return tuple([a]+beta+gam)
def canon(c): return min(coef_perm(c,p) for p in PERMS)
def valid_exact(c):
    vals=[c[0]+sum(c[1+i]*int(z[i]) for i in range(5))+sum(c[6+e]*int(z[i])*int(z[j]) for e,(i,j) in enumerate(EDGES)) for z in Z]
    return all(g>=0 for g in c[6:]) and all(min(vals[m],vals[31-m])>=0 and max(vals[m],vals[31-m])>=1 for m in range(32))
base=[tuple(Fr(v) for v in c.coefficients()) for c in make_certificates()]
newf='new_certs.json'
extra=[tuple(Fr(s) for s in c) for c in json.load(open(newf))] if os.path.exists(newf) else []
def catalog():
    S=set(base)
    for c in extra: S|={coef_perm(c,p) for p in PERMS}
    return sorted(S)
CAT=catalog(); COEF=np.array([[float(v) for v in c] for c in CAT])
print('catalog size',len(CAT),'extra orbits',len(extra))
seed=int(sys.argv[1]); secs=float(sys.argv[2]); rng=np.random.default_rng(seed); T0=time.time()
inR=lambda q: q.max()<0.5 and q.min()>0 and q.sum()-3*q.min()<=1+1e-12
n=0; hits=0; newfound=[]; examples=[]
while time.time()-T0<secs:
    k=rng.integers(2,12); S=rng.choice(32,k,replace=False); w=rng.dirichlet(np.ones(k)*rng.choice([0.3,1,3]))
    # rational pi with denominator 2000
    cnt=np.floor(w*2000).astype(int); cnt[0]+=2000-cnt.sum()
    if (cnt<0).any(): continue
    piF=[Fr(0)]*32
    for s,cn in zip(S,cnt): piF[s]=Fr(int(cn),2000)
    pi=np.array([float(v) for v in piF]); q=Z.T@pi
    if not inR(q) or q.sum()<=1: continue
    J=np.array([pi@(Z[:,i]*Z[:,j]) for i,j in EDGES]); cap=np.array([min(q[i],q[j]) for i,j in EDGES])
    t=np.round((rng.uniform(0,1,P)**3*(rng.uniform(size=P)<0.5))*1000)/1000
    if rng.integers(3)==0: t[:]=0
    b=J+(cap-J)*t; n+=1
    u=(COEF@np.r_[1,q,b]).min(); r,x=rstar(q,b)
    if u>r+1e-6:
        c=tuple(Fr(round(v*6))/6 for v in x)   # snap to sixths
        if not valid_exact(c): print('snap failed',np.round(x,4)); continue
        # exact recheck at rational theta
        qF=[sum(piF[m]*int(Z[m,i]) for m in range(32)) for i in range(5)]
        JF=[sum(piF[m]*int(Z[m,i])*int(Z[m,j]) for m in range(32)) for i,j in EDGES]
        bF=[JF[e]+(min(qF[i],qF[j])-JF[e])*Fr(str(t[e])) for e,(i,j) in enumerate(EDGES)]
        assert all(0<v<Fr(1,2) for v in qF) and sum(qF)-3*min(qF)<=1
        th=[Fr(1)]+qF+bF
        LQ=sum(a*b_ for a,b_ in zip(c,th)); UF=min(sum(a*b_ for a,b_ in zip(cc,th)) for cc in CAT)
        if LQ<UF:
            hits+=1; cn=canon(c)
            if cn not in newfound: newfound.append(cn); examples.append(dict(cert=[str(v) for v in cn],gap=float(UF-LQ),q=[str(v) for v in qF],b=[str(v) for v in bF],U=str(UF),LQ=str(LQ)))
print(f'R\\R0 instances {n}; exact-verified improvements {hits}; new orbits {len(newfound)}; time {time.time()-T0:.0f}s')
for e in examples: print(e['cert'],'gap',round(e['gap'],4))
json.dump(examples,open(f'examples_{seed}_{len(extra)}.json','w'))
if len(sys.argv)>3:  # append to catalog
    allx=[tuple(str(v) for v in c) for c in extra]+[tuple(str(v) for v in c) for c in newfound]
    json.dump(allx,open(newf,'w')); print('saved',len(allx),'extra orbits')
