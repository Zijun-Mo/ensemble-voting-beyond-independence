"""Walsh/Chow decomposition of the accuracy of a self-dual rule:
acc(g) = 1/2 + 1/2 sum_i Ghat(i)(2p_i-1) + 1/2 sum_{|S|>=3 odd} Ghat(S) E[chi_S(s)],  s_i = 2Z_i-1, G = 2g-1.
Degree-1 part identified by p; even-degree Walsh moments identified by unlabeled data; odd degree>=3 not identified.
Also: dimension of E(m,p) = 2^{K-1}-K and of P(p,r) = 2^K-1-K-C(K,2); level-1 Walsh weight of rules."""
import itertools, json
import numpy as np
from core import *
rng = np.random.default_rng(5)
out = {}
for K in (3, 5):
    Z = patterns(K); S = 2 * Z - 1
    subsets = [c for k in range(K + 1) for c in itertools.combinations(range(K), k)]
    chi = np.array([np.prod(S[:, list(c)], axis=1) if c else np.ones(2 ** K) for c in subsets]).T  # 2^K x 2^K
    G, W = enumerate_selfdual_ltfs(K, n_random=20000)
    err = 0
    for _ in range(200):
        pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3))
        p = marginals(pi, K)
        g = G[rng.integers(len(G))].astype(float); Gs = 2 * g - 1
        ghat = chi.T @ Gs / 2 ** K
        mom = chi.T @ pi
        deg = np.array([len(c) for c in subsets])
        assert np.all(np.abs(ghat[deg % 2 == 0]) < 1e-12)
        acc = 0.5 + 0.5 * sum(ghat[k] * (2 * p[subsets[k][0]] - 1) for k in range(len(subsets)) if deg[k] == 1) \
              + 0.5 * sum(ghat[k] * mom[k] for k in range(len(subsets)) if deg[k] >= 3 and deg[k] % 2 == 1)
        err = max(err, abs(acc - g @ pi))
    out[f"K{K}_walsh_identity_maxerr"] = err
    # level-1 weight of each rule type
    lw = {}
    for k in range(len(G)):
        Gs = 2 * G[k] - 1; ghat = chi.T @ Gs / 2 ** K
        deg = np.array([len(c) for c in subsets])
        w = tuple(sorted(canonical_weights(G[k], K, maxw=5), reverse=True))
        lw[str(w)] = round(float(np.sum(ghat[deg == 1] ** 2)), 4)
    out[f"K{K}_level1_weight_by_type"] = lw
    # dimensions
    A_e, b_e = set_even(np.full(K, .7), np.full(2 ** K, 2 ** -(K - 1)), K)
    out[f"K{K}_dim_E"] = int(2 ** K - np.linalg.matrix_rank(A_e))
    A_p, b_p = set_pairwise(np.full(K, .7), {(i, j): .5 for i in range(K) for j in range(i + 1, K)}, K)
    out[f"K{K}_dim_P"] = int(2 ** K - np.linalg.matrix_rank(A_p))
ps = np.array([0.55, 0.6, 2/3, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95])
out["K3_rho_star"] = {f"{p:.3f}": round(float(((5 * p - 2) / 3 - p * p) / (p * (1 - p))), 4) for p in ps}
print(json.dumps(out, indent=1))
json.dump(out, open("results/q1b_walsh.json", "w"), indent=2)
