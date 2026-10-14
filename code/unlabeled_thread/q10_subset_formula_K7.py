"""Q10 (completion session).
(1) K=5, pairwise information: check   maxmin over all 81 self-dual LTFs
        == max{ max_i p_i ,  max_{triples T} A3(T) [K=3 closed form on the triple's own (p,r)],  L_5(majority) }.
(2) K=7, pairwise information: stress test of 'subset majorities suffice' -- compare the best odd-subset majority
    (7+35+21+1 = 64 rules) with (a) 400 random positive weight vectors and (b) local search over weights
    maximising the worst-case accuracy, on random + sparse-Dirichlet + anti laws.
(3) K=7 with full unlabeled law E(m,p): same comparison (expect weighted rules can win).
"""
import json, itertools
import numpy as np
from core import *
from q2b_regret_and_anticorr import anti_instance

rng = np.random.default_rng(777)
out = {}


def A3(p3, r3):
    e = 1 - p3
    T = {k: 1 - p3[k[0]] - p3[k[1]] + r3[k] for k in r3}
    D = {k: e[k[0]] + e[k[1]] - 2 * T[k] for k in r3}
    d = [(D[(0, 1)] + D[(0, 2)] - D[(1, 2)]) / 2, (D[(0, 1)] + D[(1, 2)] - D[(0, 2)]) / 2, (D[(0, 2)] + D[(1, 2)] - D[(0, 1)]) / 2]
    return 1 - min(sum(T.values()), min(e[i] + d[i] for i in range(3)))


def gen(K, rng):
    t = rng.integers(5); Z = patterns(K)
    if t == 0:
        return anti_instance(K, rng)
    if t == 1:
        pi = rng.dirichlet(np.ones(2 ** K) * 10 ** rng.uniform(-1.5, 0.3))
        q = rng.uniform(0.55, 0.95, size=K)
        return 0.5 * pi + 0.5 * np.prod(np.where(Z == 1, q, 1 - q), axis=1)
    return family_instance(K, ["independent", "mixture", "bloc"][t - 2], rng)


# (1)
K = 5; Z = patterns(K)
G, _ = enumerate_selfdual_ltfs(K, n_random=20000); Gf = G.astype(float)
maj = (Z.sum(1) >= 3).astype(float)
errs, n = [], 0
while n < 300:
    pi = gen(K, rng); p = marginals(pi, K)
    if np.any(p <= 0.5):
        continue
    r = pairwise_joint(pi, K); A, b = set_pairwise(p, r, K)
    mm = max(lp_min(g, A, b) for g in Gf)
    tri = []
    for T in itertools.combinations(range(K), 3):
        tri.append(A3(p[list(T)], {(a, c): r[(T[a], T[c])] for a in range(3) for c in range(a + 1, 3)}))
    formula = max(p.max(), max(tri), lp_min(maj, A, b))
    errs.append(abs(mm - formula)); n += 1
out["K5_pairwise_maxmin_subset_formula"] = dict(n=n, max_abs_err=float(max(errs)))
print(out, flush=True)

# (2),(3) K=7
K = 7; Z = patterns(K); S = 2 * Z - 1
subs = [c for k in (1, 3, 5, 7) for c in itertools.combinations(range(K), k)]
subG = np.array([(Z[:, list(c)].sum(1) > len(c) / 2).astype(float) for c in subs])


def worst(w, A, b):
    m = S @ w
    if np.any(np.abs(m) < 1e-12):
        return -1
    return lp_min((m > 0).astype(float), A, b)


for info in ("pw", "full"):
    gains, n = [], 0
    while n < (40 if info == "pw" else 25):
        pi = gen(K, rng); p = marginals(pi, K)
        if np.any(p <= 0.5):
            continue
        if info == "pw":
            A, b = set_pairwise(p, pairwise_joint(pi, K), K)
        else:
            A, b = set_even(p, pi + pi[flip_index(K)], K)
        best_sub = max(lp_min(g, A, b) for g in subG)
        best_w = -1
        for _ in range(150):
            best_w = max(best_w, worst(rng.exponential(size=K) * (rng.random(K) < 0.8) + 1e-6 * rng.random(K), A, b))
        # local search from log-odds weights and from random starts
        for s in range(3):
            w = np.log(p / (1 - p)) if s == 0 else rng.exponential(size=K)
            cur = worst(w, A, b)
            for it in range(60):
                cand = np.clip(w * np.exp(0.4 * rng.normal(size=K)), 1e-6, None)
                v = worst(cand, A, b)
                if v >= cur:
                    w, cur = cand, v
            best_w = max(best_w, cur)
        gains.append(best_w - best_sub); n += 1
    g = np.array(gains)
    out[f"K7_{info}"] = dict(n=n, frac_weighted_beats_subset=float(np.mean(g > 1e-7)), max_gain=float(g.max()), mean_gain=float(g.mean()))
    print(info, out[f"K7_{info}"], flush=True)
json.dump(out, open("results/q10_results.json", "w"), indent=1)
