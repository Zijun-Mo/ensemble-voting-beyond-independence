"""Core utilities: patterns, self-dual LTF enumeration, instance generators, and LPs.

Conventions
-----------
Z in {0,1}^K : Z_i = 1 iff classifier i correct.  Y in {0,1}.  f_i = Y if Z_i=1 else 1-Y.
Patterns are indexed by integers idx in [0, 2^K); Z[idx] = bits of idx (bit i -> classifier i).
flip(idx) = idx XOR (2^K - 1)  (antipodal pattern 1-z).
alpha(z) = P(Y=1, Z=z), beta(z) = P(Y=0, Z=z), pi = alpha + beta (law of Z),
phi(v)   = P(f=v) = alpha(v) + beta(1-v)  (law of the prediction vector, from UNLABELED data).
A rule g: {0,1}^K -> {0,1} is applied to the prediction vector f.  It is self-dual iff g(1-v)=1-g(v);
for self-dual g, 1{g(f)=Y} = g(Z) pointwise, hence acc(g) = sum_z pi(z) g(z).
"""
import itertools
import numpy as np
from scipy.optimize import linprog

LP_OPTS = dict(method="highs")


def patterns(K):
    idx = np.arange(2 ** K)
    Z = ((idx[:, None] >> np.arange(K)[None, :]) & 1).astype(float)
    return Z


def flip_index(K):
    return np.arange(2 ** K) ^ (2 ** K - 1)


def reps(K):
    """Representatives of antipodal pairs: patterns with bit K-1 == 0."""
    return np.arange(2 ** (K - 1))


# ----------------------------------------------------------------------------------------
# Self-dual weighted threshold rules (positive weights, no ties)
# ----------------------------------------------------------------------------------------
def ltf_table(w, K):
    Z = patterns(K)
    S = 2 * Z - 1
    m = S @ np.asarray(w, float)
    if np.any(np.abs(m) < 1e-12):
        return None
    return (m > 0).astype(int)


def enumerate_selfdual_ltfs(K, n_random=200000, seed=0):
    """All self-dual LTFs with nonnegative weights (and no ties), by integer weights + random weights."""
    rng = np.random.default_rng(seed)
    tables = {}
    # integer weights up to 7 (covers all types for K<=5)
    for w in itertools.product(range(0, 8), repeat=K):
        if sum(w) == 0:
            continue
        t = ltf_table(w, K)
        if t is not None:
            tables.setdefault(tuple(t), np.array(w, float))
    for _ in range(n_random):
        w = rng.exponential(size=K) * (rng.random(K) < 0.85)
        if w.sum() == 0:
            continue
        t = ltf_table(w, K)
        if t is not None:
            tables.setdefault(tuple(t), w)
    keys = sorted(tables.keys())
    G = np.array(keys, int)
    W = np.array([tables[k] for k in keys])
    return G, W


def canonical_weights(table, K, maxw=9):
    """Smallest integer weight vector realising a given truth table (for display)."""
    best = None
    for w in itertools.product(range(0, maxw + 1), repeat=K):
        if sum(w) == 0 or sum(w) % 2 == 0:
            continue
        t = ltf_table(w, K)
        if t is not None and np.array_equal(t, table):
            if best is None or sum(w) < sum(best):
                best = w
    return best


# ----------------------------------------------------------------------------------------
# Instance generators
# ----------------------------------------------------------------------------------------
def mixture_of_products(K, rng, n_comp=3, acc_range=(0.6, 0.85), diff_scale=1.5, noise=0.6):
    """Latent-class (mixture of product Bernoulli) law of Z with shared 'difficulty' factor.
    Returns exact pi over {0,1}^K."""
    Z = patterns(K)
    a = np.log(np.array(rng.uniform(*acc_range, size=K)) / (1 - np.array(rng.uniform(*acc_range, size=K))))
    a = np.log(rng.uniform(*acc_range, size=K) / (1 - rng.uniform(*acc_range, size=K)))
    b = rng.uniform(0.0, diff_scale, size=K)
    wts = rng.dirichlet(np.ones(n_comp))
    pi = np.zeros(2 ** K)
    for m in range(n_comp):
        u = rng.normal()
        q = 1 / (1 + np.exp(-(a + b * u + noise * rng.normal(size=K))))
        pi += wts[m] * np.prod(np.where(Z == 1, q, 1 - q), axis=1)
    return pi / pi.sum()


def random_joint(K, rng, eta=None, **kw):
    """Random joint law of (Y,Z): class-conditional laws differ.  Returns alpha, beta, eta."""
    if eta is None:
        eta = rng.uniform(0.25, 0.75)
    p1 = mixture_of_products(K, rng, **kw)
    p0 = mixture_of_products(K, rng, **kw)
    return eta * p1, (1 - eta) * p0, eta


def phi_from(alpha, beta, K):
    fl = flip_index(K)
    return alpha + beta[fl]


def marginals(pi, K):
    return patterns(K).T @ pi


def pairwise_joint(pi, K):
    Z = patterns(K)
    return {(i, j): float(pi @ (Z[:, i] * Z[:, j])) for i in range(K) for j in range(i + 1, K)}


def pair_masses(phi, K):
    """m(z) = phi(z)+phi(1-z) = pi(z)+pi(1-z) (identified 'even part' times 2)."""
    return phi + phi[flip_index(K)]


# ----------------------------------------------------------------------------------------
# Identified sets as {x >= 0 : A x = b} (+ optional inequalities)
# ----------------------------------------------------------------------------------------
def set_pairwise(p, r, K):
    """pi >= 0, sum pi = 1, marginals p, pairwise joint successes r[(i,j)]."""
    Z = patterns(K)
    rows, rhs = [np.ones(2 ** K)], [1.0]
    for i in range(K):
        rows.append(Z[:, i]); rhs.append(p[i])
    for i in range(K):
        for j in range(i + 1, K):
            rows.append(Z[:, i] * Z[:, j]); rhs.append(r[(i, j)])
    return np.array(rows), np.array(rhs)


def set_even(p, m, K):
    """pi >= 0, pi(z)+pi(1-z) = m(z) for representatives, marginals p.  (Case (a), prior unknown.)"""
    Z = patterns(K)
    fl = flip_index(K)
    rows, rhs = [], []
    for v in reps(K):
        row = np.zeros(2 ** K); row[v] += 1; row[fl[v]] += 1
        rows.append(row); rhs.append(m[v])
    for i in range(K):
        rows.append(Z[:, i]); rhs.append(p[i])
    return np.array(rows), np.array(rhs)


def lp_range(c, A, b, A_ub=None, b_ub=None, nvar=None, bounds=None):
    """min and max of c.x over {x>=0, Ax=b, A_ub x <= b_ub}."""
    out = []
    for sgn in (1, -1):
        res = linprog(sgn * c, A_eq=A, b_eq=b, A_ub=A_ub, b_ub=b_ub,
                      bounds=bounds if bounds is not None else (0, None), **LP_OPTS)
        if res.status != 0:
            return None
        out.append(sgn * res.fun)
    return out[0], out[1]


def lp_min(c, A, b, bounds=(0, None), A_ub=None, b_ub=None):
    res = linprog(c, A_eq=A, b_eq=b, A_ub=A_ub, b_ub=b_ub, bounds=bounds, **LP_OPTS)
    if res.status != 0:
        return None
    return res.fun


# ----------------------------------------------------------------------------------------
# Direct model-level LP over (alpha, beta): the ground truth for identification questions.
# ----------------------------------------------------------------------------------------
def direct_ab_range(cz, phi, p, K, eta=None, cz_beta=None):
    """Range of sum_z alpha(z) cz(z) + beta(z) cz_beta(z) over all joint laws (alpha,beta) >= 0
    consistent with phi (law of f), marginals p (of pi = alpha+beta), and optionally P(Y=1)=eta.
    For self-dual rules use cz_beta = cz (objective = sum pi g)."""
    if cz_beta is None:
        cz_beta = cz
    n = 2 ** K
    Z = patterns(K)
    fl = flip_index(K)
    rows, rhs = [], []
    for v in range(n):  # alpha(v) + beta(1-v) = phi(v)
        row = np.zeros(2 * n); row[v] = 1; row[n + fl[v]] = 1
        rows.append(row); rhs.append(phi[v])
    for i in range(K):
        row = np.concatenate([Z[:, i], Z[:, i]]); rows.append(row); rhs.append(p[i])
    if eta is not None:
        row = np.concatenate([np.ones(n), np.zeros(n)]); rows.append(row); rhs.append(eta)
    A, b = np.array(rows), np.array(rhs)
    c = np.concatenate([cz, cz_beta])
    return lp_range(c, A, b)


def tv_range(cz, phi, p, K, eta):
    """Range of sum pi g over the even-part set + TV(pi,phi) <= 1-eta + TV(pi, phi o flip) <= eta.
    Variables: pi (n), u (n) >= |pi - phi|, w (n) >= |pi - phi~|."""
    n = 2 ** K
    fl = flip_index(K)
    A_e, b_e = set_even(p, pair_masses(phi, K), K)
    A_eq = np.hstack([A_e, np.zeros((A_e.shape[0], 2 * n))])
    I = np.eye(n)
    Zn = np.zeros((n, n))
    phit = phi[fl]
    A_ub = np.vstack([
        np.hstack([I, -I, Zn]), np.hstack([-I, -I, Zn]),        # |pi - phi| <= u
        np.hstack([I, Zn, -I]), np.hstack([-I, Zn, -I]),        # |pi - phit| <= w
        np.concatenate([np.zeros(n), 0.5 * np.ones(n), np.zeros(n)])[None, :],
        np.concatenate([np.zeros(n), np.zeros(n), 0.5 * np.ones(n)])[None, :],
    ])
    b_ub = np.concatenate([phi, -phi, phit, -phit, [1 - eta], [eta]])
    c = np.concatenate([cz, np.zeros(2 * n)])
    return lp_range(c, A_eq, b_e, A_ub=A_ub, b_ub=b_ub)


# ----------------------------------------------------------------------------------------
# Max-min over mixtures of rules / over all randomized predictors
# ----------------------------------------------------------------------------------------
def maxmin_mixture(G, A, b):
    """max_{omega in simplex} min_{x>=0, Ax=b} sum_k omega_k G_k . x   (G: n_rules x n)."""
    nr, n = G.shape
    nd = A.shape[0]
    # variables: y (nd, free), omega (nr, >=0); maximize b.y  s.t. A^T y - G^T omega <= 0, sum omega = 1
    c = np.concatenate([-b, np.zeros(nr)])
    A_ub = np.hstack([A.T, -G.T])
    b_ub = np.zeros(n)
    A_eq = np.concatenate([np.zeros(nd), np.ones(nr)])[None, :]
    bounds = [(None, None)] * nd + [(0, None)] * nr
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=bounds, **LP_OPTS)
    if res.status != 0:
        return None, None
    return -res.fun, res.x[nd:]


def maxmin_selfdual_randomized(A, b, K):
    """max over randomized self-dual g (g(z)+g(1-z)=1, g in [0,1]) of min_{x in set} g.x."""
    n = 2 ** K
    fl = flip_index(K)
    nd = A.shape[0]
    # variables y (nd free), g (n in [0,1]); max b.y s.t. A^T y - g <= 0; g(z)+g(flip z)=1
    c = np.concatenate([-b, np.zeros(n)])
    A_ub = np.hstack([A.T, -np.eye(n)])
    rows = []
    for v in reps(K):
        row = np.zeros(nd + n); row[nd + v] = 1; row[nd + fl[v]] = 1; rows.append(row)
    A_eq = np.array(rows)
    bounds = [(None, None)] * nd + [(0, 1)] * n
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(n), A_eq=A_eq, b_eq=np.ones(len(rows)), bounds=bounds, **LP_OPTS)
    if res.status != 0:
        return None, None
    return -res.fun, res.x[nd:]


def bf_value(phi, p, K, eta=None):
    """Balsubramani-Freund-type population game with EQUALITY accuracy constraints:
    max over randomized predictors g:{0,1}^K->[0,1] (any, not nec. self-dual) of the worst-case
    accuracy over q_v = P(Y=1|f=v) in [0,1] consistent with p (and eta)."""
    n = 2 ** K
    Z = patterns(K)
    S = 2 * Z - 1
    d = p - phi @ (1 - Z)  # sum_v y_v s_i(v) = d_i
    rows_extra = []
    # variables: g (n, [0,1]), lam (K free), [nu free if eta], mu (n, >=0)
    ne = 1 if eta is not None else 0
    nv = n + K + ne + n
    # objective: maximize sum phi(1-g) + lam.d + nu*eta - phi.mu
    c = np.zeros(nv)
    c[:n] = phi  # minimize -(...) -> +phi.g  (constant sum phi dropped)
    c[n:n + K] = -d
    if ne:
        c[n + K] = -eta
    c[n + K + ne:] = phi
    # constraint: s(v).lam + nu - mu_v - 2 g_v <= -1
    A_ub = np.zeros((n, nv))
    A_ub[:, :n] = -2 * np.eye(n)
    A_ub[:, n:n + K] = S
    if ne:
        A_ub[:, n + K] = 1.0
    A_ub[:, n + K + ne:] = -np.eye(n)
    b_ub = -np.ones(n)
    bounds = [(0, 1)] * n + [(None, None)] * (K + ne) + [(0, None)] * n
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, **LP_OPTS)
    if res.status != 0:
        return None, None
    val = phi.sum() - res.fun
    return val, res.x[:n]


# ----------------------------------------------------------------------------------------
# Dual ("unlabeled-margin") formula for the worst-case accuracy of a self-dual rule
# ----------------------------------------------------------------------------------------
def dual_value(lam, g, phi, p, K, upper=False):
    """(1/2)[1 + lam.(2p-1) -/+ E_phi |1 - G(f) lam.s(f)|],  G = 2g-1 in {+-1}."""
    Z = patterns(K)
    S = 2 * Z - 1
    Gs = 2 * g - 1
    term = phi @ np.abs(1 - Gs * (S @ lam))
    base = 1 + lam @ (2 * p - 1)
    return 0.5 * (base + term) if upper else 0.5 * (base - term)


# ---------------------------------------------------------------------------
# Instance families used in Q2 (copied from q2_robust_weights.py)
# ---------------------------------------------------------------------------
def family_instance(K, fam, rng):
    Z = patterns(K)
    if fam == "independent":
        p = rng.uniform(0.55, 0.9, size=K)
        pi = np.prod(np.where(Z == 1, p, 1 - p), axis=1)
    elif fam == "mixture":
        pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3), acc_range=(0.55, 0.9))
    elif fam == "bloc":
        # first nb classifiers share a latent correct/incorrect draw with prob c; others independent
        nb = 2 if K == 3 else 3
        c = rng.uniform(0.3, 0.95)
        pb = rng.uniform(0.6, 0.85)
        po = rng.uniform(0.55, 0.85, size=K - nb)
        pi = np.zeros(2 ** K)
        for zb in (0, 1):  # common draw
            pc = pb if zb == 1 else 1 - pb
            q = np.concatenate([np.full(nb, c * zb + (1 - c) * pb), po])
            pi += pc * np.prod(np.where(Z == 1, q, 1 - q), axis=1)
    pi = np.clip(pi, 0, None)
    return pi / pi.sum()


