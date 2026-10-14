"""Q6 (completion session).
(1) Adversarial search for a counterexample to the conjecture
      'max-min over self-dual weighted threshold rules (given p + pairwise, or p + full unlabeled law)
       is always attained by an UNWEIGHTED odd-subset majority'.
    K=4: the only non-subset type is (2,1,1,1).  K=5: types (2,1,1,1,0), (3,1,1,1,1), (2,2,1,1,1), (3,2,2,1,1).
    Search: random Dirichlet laws (small concentration -> sparse / extreme dependence) + hill climbing on pi
    to maximise  gain = maxmin(all rules) - maxmin(subset majorities).
(2) K=5 majority: harvest optimal dual vertices of the sharp-lower-bound LP over P(p,r) (pairwise info),
    dedupe up to S_5 symmetry, and test whether max over the harvested certificate family reproduces the
    LP on fresh instances (=> numerically verified closed form 'max of N Bonferroni-type inequalities').
"""
import json, itertools
import numpy as np
from scipy.optimize import linprog
from core import *
from q2b_regret_and_anticorr import anti_instance

rng = np.random.default_rng(31337)
out = {}


def rules(K):
    G, W = enumerate_selfdual_ltfs(K, n_random=20000)
    is_subset = []
    for k in range(len(G)):
        w = canonical_weights(G[k], K)
        is_subset.append(w is not None and set(w) <= {0, 1})
    return G.astype(float), np.array(is_subset)


def gain(pi, K, Gf, sub, info):
    p = marginals(pi, K)
    if info == "pw":
        A, b = set_pairwise(p, pairwise_joint(pi, K), K)
    else:
        A, b = set_even(p, pi + pi[flip_index(K)], K)
    lows = np.array([lp_min(g, A, b) for g in Gf])
    return lows.max() - lows[sub].max(), lows


for K in (4, 5):
    Gf, sub = rules(K)
    out[f"K{K}_n_rules"] = int(len(Gf)); out[f"K{K}_n_subset"] = int(sub.sum())
    for info in ("pw", "full"):
        best_gain, n_eval = -1, 0
        # random restarts
        n_starts = 60 if K == 4 else 25
        for s in range(n_starts):
            conc = 10 ** rng.uniform(-1.5, 0.5)
            pi = rng.dirichlet(np.ones(2 ** K) * conc)
            # bias toward p>1/2 by mixing with a good product law
            q = rng.uniform(0.55, 0.95, size=K); Zk = patterns(K)
            pi = 0.5 * pi + 0.5 * np.prod(np.where(Zk == 1, q, 1 - q), axis=1)
            gcur, _ = gain(pi, K, Gf, sub, info); n_eval += 1
            step = 0.3
            for it in range(40 if K == 4 else 20):
                cand = pi * np.exp(step * rng.normal(size=2 ** K)); cand /= cand.sum()
                if np.rng if False else False:
                    pass
                gc, _ = gain(cand, K, Gf, sub, info); n_eval += 1
                if gc >= gcur:
                    pi, gcur = cand, gc
                step *= 0.97
            best_gain = max(best_gain, gcur)
        out[f"K{K}_{info}_adversarial"] = dict(n_lp_instances=n_eval, max_gain=float(best_gain))
        print(K, info, out[f"K{K}_{info}_adversarial"], flush=True)

# ---------------------------------------------------------------- (2) K=5 majority certificates
K = 5
Zk = patterns(K); maj = (Zk.sum(1) >= 3).astype(float)
pairs = [(i, j) for i in range(K) for j in range(i + 1, K)]
perms = list(itertools.permutations(range(K)))


def dual_of_min(c, A, b):
    res = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    return res.fun, res.eqlin.marginals  # y with  A^T y <= c, value = b.y


def canon(y):
    """canonical form under S_5 acting on (y0, y_i, Y_ij)."""
    best = None
    for pm in perms:
        inv = np.argsort(pm)
        yi = y[1:1 + K][list(pm)]
        Y = np.zeros((K, K))
        for k, (i, j) in enumerate(pairs):
            Y[i, j] = Y[j, i] = y[1 + K + k]
        Yp = Y[np.ix_(pm, pm)]
        vec = np.concatenate([[y[0]], yi, [Yp[i, j] for (i, j) in pairs]])
        key = tuple(np.round(vec, 6))
        if best is None or key < best:
            best = key
    return best


def gen(rng):
    t = rng.integers(5)
    if t == 0:
        pi = anti_instance(K, rng)
    elif t == 1:
        pi = rng.dirichlet(np.ones(32) * 10 ** rng.uniform(-1.5, 0.5))
        q = rng.uniform(0.55, 0.95, size=K)
        pi = 0.5 * pi + 0.5 * np.prod(np.where(Zk == 1, q, 1 - q), axis=1)
    else:
        pi = family_instance(K, ["independent", "mixture", "bloc"][t - 2], rng)
    return pi


certs = {}
for it in range(1500):
    pi = gen(rng)
    p = marginals(pi, K)
    A, b = set_pairwise(p, pairwise_joint(pi, K), K)
    val, y = dual_of_min(maj, A, b)
    if np.max(A.T @ y - maj) > 1e-7:
        continue
    certs.setdefault(canon(y), y)
out["K5_maj_n_distinct_cert_orbits_harvested"] = len(certs)
# expand orbits
Ys = []
for key in certs:
    y = np.array(key)
    for pm in perms:
        yi = y[1:1 + K]; Y = np.zeros((K, K))
        for k, (i, j) in enumerate(pairs):
            Y[i, j] = Y[j, i] = y[1 + K + k]
        newyi = np.zeros(K); newY = np.zeros((K, K))
        newyi[list(pm)] = yi
        newY[np.ix_(pm, pm)] = Y
        Ys.append(np.concatenate([[y[0]], newyi, [newY[i, j] for (i, j) in pairs]]))
Ys = np.unique(np.round(np.array(Ys), 9), axis=0)
out["K5_maj_n_certificates_total"] = int(len(Ys))
# test on fresh instances
errs = []
for it in range(1000):
    pi = gen(rng)
    p = marginals(pi, K)
    A, b = set_pairwise(p, pairwise_joint(pi, K), K)
    lp = lp_min(maj, A, b)
    closed = float(np.max(Ys @ b))
    errs.append(lp - closed)
errs = np.array(errs)
out["K5_maj_certfamily_fresh"] = dict(n=1000, frac_exact=float(np.mean(np.abs(errs) < 1e-7)), max_gap=float(errs.max()),
                                      min_gap=float(errs.min()))
print(out["K5_maj_certfamily_fresh"], len(certs), len(Ys))
# describe orbit representatives in 'tandem' coordinates: bound = y0 + sum y_i p_i + sum Y_ij r_ij
reps_list = []
for key in certs:
    y = np.array(key)
    reps_list.append(dict(y0=y[0], y_i=list(y[1:1 + K]), Y_ij=list(y[1 + K:])))
out["K5_maj_cert_orbit_representatives"] = reps_list
json.dump(out, open("results/q6_results.json", "w"), indent=1, default=float)
