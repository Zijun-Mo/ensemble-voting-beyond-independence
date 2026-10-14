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


def _exact_solve(rows, rhs, ncols):
    """Solve rows * x = rhs exactly (Fractions, Gauss-Jordan); return one solution or None if inconsistent."""
    M = [list(r) + [v] for r, v in zip(rows, rhs)]
    piv = []; r0 = 0
    for c in range(ncols):
        p = next((r for r in range(r0, len(M)) if M[r][c] != 0), None)
        if p is None: continue
        M[r0], M[p] = M[p], M[r0]
        inv = 1 / M[r0][c]; M[r0] = [x * inv for x in M[r0]]
        for r in range(len(M)):
            if r != r0 and M[r][c] != 0:
                f = M[r][c]; M[r] = [a - f * b_ for a, b_ in zip(M[r], M[r0])]
        piv.append(c); r0 += 1
    if any(all(x == 0 for x in row[:-1]) and row[-1] != 0 for row in M): return None
    x = [F(0)] * ncols
    for k, c in enumerate(piv): x[c] = M[k][-1]
    return x


def certified_lower_bound(K, q_minus, q_plus, b):
    """Exact rational LOWER bound on C_full (added in review, 2026-10-14).

    For every distribution pi in the empirical ambiguity set P_hat and every deterministic self-dual F,
    E_pi F >= sum_z min(pi_z, pi_{1-z}) (z over complementary pairs). Hence for ANY exactly feasible pi,
        C_full = min_F max_pi E_pi F >= max_pi min_F E_pi F >= L(pi) := sum_z min(pi_z, pi_{1-z}).
    The best such bound is the randomized value, which by LP duality equals the LP relaxation (t_z in [0,1]) of the
    no-big-M MILP. We solve that LP with dual simplex, rebuild the optimal vertex EXACTLY from its active constraints
    and positive variables (Fraction Gauss-Jordan), and re-check feasibility in Fractions. If reconstruction fails,
    a weaker but exact bound is returned (never an unverified number).
    """
    from scipy.optimize import linprog
    ql = tuple(map(fraction, q_minus)); qu = tuple(map(fraction, q_plus)); be = tuple(map(fraction, b))
    edges, reps = _items(K); N = 1 << K; allmask = N - 1; t = len(reps); nv = N + t
    X = [[(m >> i) & 1 for i in range(K)] for m in range(N)]
    rows = []; rhs = []
    for i in range(K):
        rows.append([X[m][i] for m in range(N)] + [0] * t); rhs.append(qu[i])
        rows.append([-X[m][i] for m in range(N)] + [0] * t); rhs.append(-ql[i])
    for (i, j) in edges:
        rows.append([X[m][i] * X[m][j] for m in range(N)] + [0] * t); rhs.append(be[edges.index((i, j))])
    for k, z in enumerate(reps):
        for m in (z, allmask ^ z):
            r = [0] * nv; r[N + k] = 1; r[m] = -1; rows.append(r); rhs.append(F(0))
    eq = [1] * N + [0] * t
    res = linprog(np.r_[np.zeros(N), -np.ones(t)], A_ub=np.array(rows, float), b_ub=np.array([float(v) for v in rhs]),
                  A_eq=np.array([eq], float), b_eq=[1.0], bounds=[(0, None)] * nv, method='highs-ds')
    def check(x):
        if x is None or any(v < 0 for v in x) or sum(x[:N]) != 1: return False
        return all(sum(a * v for a, v in zip(r, x) if a) <= h for r, h in zip(rows, rhs))
    if res.status == 0:
        xf = res.x; slack = np.array([float(h) for h in rhs]) - np.array(rows, float) @ xf
        act = [k for k in range(len(rows)) if abs(slack[k]) < 1e-9]
        pos = [c for c in range(nv) if xf[c] > 1e-12]
        sub = [[F(rows[k][c]) for c in pos] for k in act] + [[F(eq[c]) for c in pos]]
        y = _exact_solve(sub, [rhs[k] for k in act] + [F(1)], len(pos))
        if y is not None:
            x = [F(0)] * nv
            for c, v in zip(pos, y): x[c] = v
            if check(x):
                pi = x[:N]
                return {'certified_lower': sum(min(pi[z], pi[allmask ^ z]) for z in reps), 'pi': pi,
                        'method': 'exact vertex reconstruction', 'lp_value_float': -res.fun}
    return {'certified_lower': F(0), 'pi': None, 'method': 'fallback (trivial bound 0)', 'lp_value_float': None}


def certified_gap(K, q_minus, q_plus, b, time_limit=30.0):
    """Verified upper certificate U (rule + exact dual) and exact lower bound L; eta := U - L bounds the
    optimization error  U - C_full  rigorously."""
    up = solve_and_repair_upper(K, q_minus, q_plus, b, time_limit=time_limit)
    lo = certified_lower_bound(K, q_minus, q_plus, b)
    return {'upper': up['verified_upper'], 'lower': lo['certified_lower'],
            'eta': up['verified_upper'] - lo['certified_lower'], 'upper_cert': up, 'lower_cert': lo}


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
    print('NOTE: global MILP optimality is not certified; the exact gap bound below replaces it.')
    g=certified_gap(K,ql,qu,be)
    print('certified lower bound =',g['lower'],' eta = U - L =',g['eta'])
    # Independent voters q=1/5, b=1/25: deterministic value 3/25, randomized value 1/10 (doc result 9(b)).
    q=[F(1,5)]*K; bb=[F(1,25)]*10
    g=certified_gap(K,q,q,bb)
    print('q=1/5, b=1/25: U =',g['upper'],' L =',g['lower'],' eta =',g['eta'],' (deterministic 3/25 vs randomized 1/10)')
