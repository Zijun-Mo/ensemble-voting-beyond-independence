"""Symmetric parameters q_i = q, b_ij = b: compare exact R*_det (no-big-M MILP over
all self-dual rules) with the best 'sub-committee' rule:
  m odd  : majority of m members (others ignored)
  m even : weights (2,1,...,1) on m members (majority with one tie-breaker)
Each rule's worst case is an exact LP over all 2^K joint error patterns.
Usage: python3 sym_explore.py K n_points seconds
"""
import sys, time, itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

K = int(sys.argv[1]); NPTS = int(sys.argv[2]); TLIM = float(sys.argv[3])
E = list(itertools.combinations(range(K), 2)); M = len(E)
Z = np.array([[(s >> i) & 1 for i in range(K)] for s in range(1 << K)])
PAIR = np.array([[z[i] * z[j] for (i, j) in E] for z in Z])
ALL = (1 << K) - 1


def worst(F, q, b):
    """max_pi E_pi F  s.t. marginals = q, pair moments <= b."""
    A_eq = np.vstack([np.ones(1 << K), Z.T]); b_eq = np.r_[1.0, np.full(K, q)]
    res = linprog(-F, A_ub=PAIR.T, b_ub=np.full(M, b), A_eq=A_eq, b_eq=b_eq,
                  bounds=(0, None), method='highs')
    return -res.fun if res.status == 0 else None


def feasible(q, b):
    return worst(np.zeros(1 << K), q, b) is not None


def rstar_det(q, b):
    reps = [s for s in range(1 << K) if not (s & 1)]
    t = len(reps); n = 1 + K + M + t
    c = np.zeros(n); c[0] = 1; c[1:1 + K] = q; c[1 + K:1 + K + M] = b
    A = np.zeros((2 * t, n)); lo = np.zeros(2 * t)
    for j, r in enumerate(reps):
        for row, s, coef in ((2 * j, r, -1), (2 * j + 1, ALL ^ r, 1)):
            A[row, 0] = 1; A[row, 1:1 + K] = Z[s]; A[row, 1 + K:1 + K + M] = PAIR[s]
            A[row, 1 + K + M + j] = coef
        lo[2 * j + 1] = 1
    lb = np.r_[-np.inf, np.full(K, -np.inf), np.zeros(M), np.zeros(t)]
    ub = np.r_[np.inf, np.full(K, np.inf), np.full(M, np.inf), np.ones(t)]
    integ = np.r_[np.zeros(1 + K + M), np.ones(t)]
    res = milp(c, integrality=integ, bounds=Bounds(lb, ub),
               constraints=LinearConstraint(A, lo, np.inf),
               options={'time_limit': 60, 'mip_rel_gap': 1e-9})
    return res.fun, res.status


def committee_rule(m):
    w = np.zeros(K); w[:m] = 1
    if m % 2 == 0:
        w[0] = 2
    tot = w.sum()
    return (Z @ w > tot / 2).astype(float), tuple(int(x) for x in w[:m])


RULES = [committee_rule(m) for m in range(1, K + 1)]
rng = np.random.default_rng(K)
T = time.time(); n = 0; worst_gap = 0; bad = []; who = {}
while n < NPTS and time.time() - T < TLIM:
    q = rng.uniform(1.0 / K * 0.7, 0.5)
    b = rng.uniform(0, q) * rng.uniform(0.3, 1.0)
    if not feasible(q, b):
        continue
    vals = [worst(F, q, b) for F, _ in RULES]
    cand = min(vals); arg = RULES[int(np.argmin(vals))][1]
    ex, st = rstar_det(q, b)
    if st != 0:
        continue
    n += 1
    gap = cand - ex; worst_gap = max(worst_gap, gap)
    who[arg] = who.get(arg, 0) + 1
    if gap > 1e-7:
        bad.append((round(q, 4), round(b, 4), round(ex, 6), round(cand, 6), arg))
print(f'K={K}: {n} feasible points in {time.time()-T:.0f}s; max(candidate - exact) = {worst_gap:.2e}')
print('  best committee rule counts:', dict(sorted(who.items(), key=lambda kv: -kv[1])))
print('  points where an asymmetric rule beats every committee rule:', len(bad))
for x in bad[:8]:
    print('   q=%s b=%s exact=%s committee=%s %s' % x)
