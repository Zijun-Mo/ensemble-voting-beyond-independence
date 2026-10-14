#!/usr/bin/env python3
"""Finite-sample minimax voting, interval-marginal MILP with rational UPPER certificates.

Requires numpy, scipy. The MILP is a numerical candidate generator ONLY.
Only the subsequent Fraction-based exact constraint checks make the returned
upper error certificate rigorous for rational input confidence endpoints.
Global MILP optimality is NOT rationally certified by this script.

Input (Python API): q_minus, q_plus each length K; b a length-m list ordered by
itertools.combinations(range(K), 2). Values may be Fraction, decimal strings,
or integers; float values are parsed through decimal str (not ideal for external
statistical endpoints, which should have been outward-rounded before calling).
"""
from __future__ import annotations
from fractions import Fraction as F
from itertools import combinations
from typing import Iterable
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint


def fraction(v):
    return v if isinstance(v, F) else F(str(v))


def _items(K):
    edges = list(combinations(range(K), 2))
    reps = [s for s in range(1 << K) if (s & 1) == 0]
    return edges, reps


def verify_rational_upper(K, q_minus, q_plus, b, decisions,
                          a, lam_plus, lam_minus, gamma, primal_pi=None):
    """Verify EVERY Q(x)>=F(x) constraint exactly, and optionally primal LP feasibility.

    decisions is dict {representative_bitmask: 0 or 1}; f(rep)=t, f(complement)=1-t.
    If primal_pi is provided, also verify this distribution is feasible for the
    chosen self-dual F and report the resulting LOWER bound on its max risk.
    Lower bound for a single F is NOT a minimax lower bound across all F.
    """
    if K < 1: raise ValueError('K must be positive')
    edges, reps = _items(K)
    ql = tuple(map(fraction, q_minus)); qu = tuple(map(fraction, q_plus))
    be = tuple(map(fraction, b))
    lp = tuple(map(fraction, lam_plus)); lm = tuple(map(fraction, lam_minus))
    ga = tuple(map(fraction, gamma)); a = fraction(a)
    assert len(ql) == len(qu) == len(lp) == len(lm) == K
    assert len(be) == len(ga) == len(edges)
    assert all(F(0) <= l <= u <= F(1) for l,u in zip(ql, qu))
    assert all(v >= 0 for v in be + lp + lm + ga)
    assert set(decisions) == set(reps) and all(decisions[r] in (0, 1) for r in reps)
    allmask = (1 << K) - 1
    def risk_at(mask):
        t = decisions.get(mask)
        if t is not None: return t
        return 1 - decisions[allmask ^ mask]
    def poly(mask):
        ans = a
        for i in range(K):
            if mask >> i & 1: ans += lp[i] - lm[i]
        for ei, (i,j) in enumerate(edges):
            if (mask >> i & 1) and (mask >> j & 1): ans += ga[ei]
        return ans
    deficits = [(mask, F(risk_at(mask)) - poly(mask)) for mask in range(1 << K)
                if poly(mask) < risk_at(mask)]
    if deficits:
        raise ValueError(f'Dual certificate infeasible; first violations {deficits[:3]}')
    upper = a + sum(x*y - w*z for x,y,w,z in zip(lp,qu,lm,ql)) + sum(x*y for x,y in zip(ga,be))
    ans = {'verified_upper':upper, 'rule_error_masks':[x for x in range(1<<K) if risk_at(x)],
           'dual_min_slack':min(poly(x) - risk_at(x) for x in range(1<<K))}
    if primal_pi is not None:
        p = tuple(map(fraction,primal_pi))
        assert len(p)==(1<<K) and all(x>=0 for x in p) and sum(p)==1
        for i in range(K):
            marginal = sum(p[x] for x in range(1<<K) if x>>i & 1)
            assert ql[i] <= marginal <= qu[i], (i,marginal,ql[i],qu[i])
        for ei,(i,j) in enumerate(edges):
            pair = sum(p[x] for x in range(1<<K) if (x>>i & 1) and (x>>j & 1))
            assert pair <= be[ei], (ei,pair,be[ei])
        ans['verified_primal_lower_for_this_rule'] = sum(p[x]*risk_at(x) for x in range(1<<K))
        assert ans['verified_primal_lower_for_this_rule'] <= upper
    return ans


def solve_and_repair_upper(K, q_minus, q_plus, b, max_denominator=10**7, time_limit=30.0):
    """MILP with no big-M; return exact-rational repaired upper-bound certificate.

    Repair: rationalize nonnegative multipliers; then increase intercept a by the
    largest exact Q>=F deficit, which guarantees feasibility without trusting
    floating-point tolerances. The objective may grow; no exact optimality claim.
    """
    ql = tuple(map(fraction,q_minus)); qu = tuple(map(fraction,q_plus))
    be = tuple(map(fraction,b))
    edges,reps=_items(K)
    m=len(edges); t=len(reps)
    # Index layout: [a, lambda_plus[K], lambda_minus[K], gamma[m], t[t]]
    n=1+2*K+m+t
    ia=0; ip=1; im=1+K; ig=1+2*K; it=ig+m
    obj=np.zeros(n)
    obj[ia]=1
    for i in range(K):obj[ip+i]=float(qu[i]);obj[im+i]=-float(ql[i])
    for e in range(m):obj[ig+e]=float(be[e])
    A=np.zeros((2*t,n)); lower=np.zeros(2*t)
    allmask=(1<<K)-1
    for j,rep in enumerate(reps):
        # Q(rep) - t_rep >= 0
        for r,mask,coef in ((2*j,rep,-1),(2*j+1,allmask^rep,1)):
            A[r,ia]=1
            for i in range(K):
                if (mask>>i)&1: A[r,ip+i]=1;A[r,im+i]=-1
            for e,(i,k) in enumerate(edges):
                if (mask>>i&1) and (mask>>k&1): A[r,ig+e]=1
            A[r,it+j]=coef
        lower[2*j+1]=1
    lb=np.zeros(n); ub=np.full(n,np.inf)
    lb[ia]=-np.inf
    for j in range(t):ub[it+j]=1
    integrality=np.zeros(n,int);integrality[it:]=1
    res=milp(obj, integrality=integrality,bounds=Bounds(lb,ub),
             constraints=LinearConstraint(A,lower,np.full(2*t,np.inf)),
             options={'time_limit':time_limit,'mip_rel_gap':0.0})
    if res.x is None:
        raise RuntimeError(f'No candidate obtained: {res.message}; check set feasibility')
    x=res.x
    decisions={reps[j]:int(round(x[it+j])) for j in range(t)}
    rn=lambda z:F(str(float(z))).limit_denominator(max_denominator)
    aa=rn(x[ia]); lp=[max(F(0),rn(x[ip+i])) for i in range(K)]
    lm=[max(F(0),rn(x[im+i])) for i in range(K)]
    ga=[max(F(0),rn(x[ig+e])) for e in range(m)]
    # Exact symbolic intercept repair.
    def Q(mask):
        ans=aa
        for i in range(K):
            if (mask>>i)&1:ans+=lp[i]-lm[i]
        for e,(i,k) in enumerate(edges):
            if (mask>>i&1) and (mask>>k&1):ans+=ga[e]
        return ans
    allmask=(1<<K)-1
    deficit=max([F(0)]+[F((decisions[mask] if mask in decisions else
                                 1-decisions[allmask^mask]))-Q(mask)
                           for mask in range(1<<K)])
    aa+=deficit
    verified=verify_rational_upper(K,ql,qu,be,decisions,aa,lp,lm,ga)
    verified.update({'a':aa,'lambda_plus':lp,'lambda_minus':lm,'gamma':ga,
                     'decisions':decisions, 'rounding_repair':deficit,
                     'solver_objective_diagnostic':res.fun,
                     'solver_lower_bound_DIAGNOSTIC_NOT_CERTIFIED':getattr(res,'mip_dual_bound',None),
                     'solver_status':res.message})
    return verified


if __name__ == '__main__':
    # Deterministic rational demo: all voters have q=0.3, b=0.05, so 5-voter
    # majority is oracle-zero. Candidate is verified EXACTLY, not by float slack.
    K=5
    ql=[F(3,10)]*K;qu=[F(3,10)]*K
    be=[F(1,20)]*10
    r=solve_and_repair_upper(K,ql,qu,be,time_limit=20)
    print('exact candidate upper =',r['verified_upper'])
    print('repaired intercept =',r['a'], 'repair =',r['rounding_repair'])
    print('MILP optimizer message:',r['solver_status'])
    print('NOTE: global MILP optimality was not independently certified.')
