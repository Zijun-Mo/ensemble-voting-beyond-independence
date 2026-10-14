"""Core utilities: robust (worst-case over dependence) accuracy of voting rules.

Z in {0,1}^K, Z_i = 1{classifier i correct}.  A rule's correctness is a Boolean
function g(Z) (for self-dual label rules, correct iff g(Z)=1).
V(g; C) = min_{law of Z in C} P(g(Z)=1), C = marginals (+ optional pairwise constraints).
"""
import itertools
import numpy as np
from scipy.optimize import linprog


def states(K):
    return np.array(list(itertools.product([0, 1], repeat=K)), dtype=np.int64)


def idx_of(Z):
    K = Z.shape[1]
    return Z @ (1 << np.arange(K - 1, -1, -1))


# ---------------------------------------------------------------- rule classes
def weighted_rule(w, Z):
    """g(z)=1 iff sum_i w_i (2 z_i - 1) > 0  (ties -> 0, avoid with generic w)."""
    return ((2 * Z - 1) @ np.asarray(w, float) > 1e-12).astype(np.int64)


def selfdual_functions(K, monotone_only=True, unanimity=True):
    """All self-dual Boolean functions on {0,1}^K (as 0/1 vectors over states(K)).
    Self-dual: g(1-z) = 1 - g(z).  Optionally keep only monotone ones."""
    Z = states(K)
    n = len(Z)
    idx = idx_of(Z)
    comp = idx_of(1 - Z)
    reps = [i for i in range(n) if Z[i, 0] == 0]  # one of each complementary pair
    nr = len(reps)
    # covers for monotonicity: z -> z + e_j
    covers = [(i, idx[i] + (1 << (K - 1 - j))) for i in range(n) for j in range(K) if Z[i, j] == 0]
    covers = np.array(covers)
    out = []
    for mask in range(1 << nr):
        g = np.zeros(n, dtype=np.int64)
        for b, r in enumerate(reps):
            v = (mask >> b) & 1
            g[r] = v
            g[comp[r]] = 1 - v
        if unanimity and g[n - 1] != 1:  # g(1...1)=1 (=> g(0..0)=0)
            continue
        if monotone_only and np.any(g[covers[:, 0]] > g[covers[:, 1]]):
            continue
        out.append(g)
    return np.array(out)


def is_threshold_selfdual(g, Z):
    """LP: exists w>=0 with sum w_i(2z_i-1) >= 1 for all winning z?  Returns weights or None."""
    K = Z.shape[1]
    Wn = Z[g == 1]
    A = -(2 * Wn - 1).astype(float)
    res = linprog(np.ones(K), A_ub=A, b_ub=-np.ones(len(Wn)), bounds=[(0, None)] * K, method="highs")
    return res.x if res.status == 0 else None


def sample_weighted_rules(K, nsamp, rng):
    """De-duplicated truth tables of weighted votes from random positive weights."""
    Z = states(K)
    seen = {}
    S = 2 * Z - 1
    gens = [
        lambda: rng.exponential(size=K),
        lambda: rng.dirichlet(0.3 * np.ones(K)),
        lambda: rng.integers(1, 8, size=K) + 1e-6 * rng.random(K),
        lambda: rng.lognormal(0, 1.5, size=K),
        lambda: np.where(rng.random(K) < 0.4, 0.0, 1.0) * rng.exponential(size=K) + 1e-9 * rng.random(K),
    ]
    for t in range(nsamp):
        w = gens[t % len(gens)]()
        s = S @ w
        if np.min(np.abs(s)) < 1e-9:
            continue
        g = (s > 0).astype(np.int64)
        key = g.tobytes()
        if key not in seen:
            seen[key] = (g, w)
    return [v[0] for v in seen.values()], [v[1] for v in seen.values()]


# ---------------------------------------------------------------- LPs
def constraint_system(p, rbar=None, req=None):
    """Feasible set C for law x of Z: sum x = 1, E Z_i = p_i,
    optional E Z_i Z_j <= rbar_ij (upper bounds) or == req_ij (exact)."""
    p = np.asarray(p, float)
    K = len(p)
    Z = states(K).astype(float)
    Aeq = [np.ones(len(Z))] + [Z[:, i] for i in range(K)]
    beq = [1.0] + list(p)
    Aub, bub = [], []
    for i in range(K):
        for j in range(i + 1, K):
            if req is not None:
                Aeq.append(Z[:, i] * Z[:, j]); beq.append(req[i][j])
            elif rbar is not None and rbar[i][j] is not None and np.isfinite(rbar[i][j]):
                Aub.append(Z[:, i] * Z[:, j]); bub.append(rbar[i][j])
    Aeq = np.array(Aeq); beq = np.array(beq)
    Aub = np.array(Aub) if Aub else None
    bub = np.array(bub) if bub else None
    return Z, Aeq, beq, Aub, bub


def V(g, p, rbar=None, req=None, return_x=False, system=None):
    """Worst-case P(g(Z)=1) over laws in C."""
    if system is None:
        system = constraint_system(p, rbar, req)
    Z, Aeq, beq, Aub, bub = system
    res = linprog(np.asarray(g, float), A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq,
                  bounds=(0, None), method="highs")
    if res.status != 0:
        return (None, None) if return_x else None
    return (res.fun, res.x) if return_x else res.fun


def robust_mixture(G, p, rbar=None, req=None, system=None):
    """max_{lambda in simplex(G)} min_{x in C} sum_g lambda_g P_x(g=1)  (one LP via duality).
    Returns (value, lambda)."""
    if system is None:
        system = constraint_system(p, rbar, req)
    Z, Aeq, beq, Aub, bub = system
    G = np.asarray(G, float)
    nG, n = G.shape
    me = Aeq.shape[0]
    mu = 0 if Aub is None else Aub.shape[0]
    # variables: lambda (nG, >=0), y (me, free), u (mu, >=0); maximise beq.y - bub.u
    c = np.concatenate([np.zeros(nG), -beq, (bub if mu else np.zeros(0))])
    # A_eq^T y - A_ub^T u - G^T lambda <= 0
    A1 = np.hstack([-G.T, Aeq.T, (-Aub.T if mu else np.zeros((n, 0)))])
    A_eq2 = np.concatenate([np.ones(nG), np.zeros(me + mu)])[None, :]
    bounds = [(0, None)] * nG + [(None, None)] * me + [(0, None)] * mu
    res = linprog(c, A_ub=A1, b_ub=np.zeros(n), A_eq=A_eq2, b_eq=[1.0], bounds=bounds, method="highs")
    if res.status != 0:
        return None, None
    return -res.fun, res.x[:nG]


def comonotone_value(g, p):
    """P(g=1) under the comonotone coupling Z_i = 1{U < p_i}."""
    p = np.asarray(p, float)
    K = len(p)
    order = np.argsort(-p, kind="stable")
    Z = states(K)
    gi = dict(zip(map(tuple, Z), g))
    val = 0.0
    # U in [p_(j+1), p_(j)) -> top-j set correct
    ps = np.concatenate([p[order], [0.0]])
    for j in range(K + 1):
        z = np.zeros(K, dtype=np.int64)
        z[order[:j]] = 1
        lo = ps[j]
        hi = 1.0 if j == 0 else ps[j - 1]
        val += gi[tuple(z)] * max(0.0, hi - lo)
    return val


def matan_lower(p):
    p = np.sort(np.asarray(p, float))[::-1]
    K = len(p); m = (K + 1) // 2
    return max([0.0] + [(p[:j].sum() - (m - 1)) / (j - m + 1) for j in range(m, K + 1)])


def rbar_from_rho(p, rho_bar):
    p = np.asarray(p, float)
    K = len(p)
    s = np.sqrt(p * (1 - p))
    R = np.full((K, K), np.nan)
    for i in range(K):
        for j in range(i + 1, K):
            R[i, j] = p[i] * p[j] + rho_bar * s[i] * s[j]
    return R
