"""Exercise certified_gap (S29_bundle) on random rational interval instances for K=5,6,7: how often the exact
lower bound reconstructs (vs. the trivial fallback 0), the certified gap eta = U - L, and timing.
Usage: python3 -I check_gap_certificate.py <S29_bundle dir>"""
import sys, time, itertools, numpy as np
from fractions import Fraction as F
sys.path.insert(0, sys.argv[1]); import task1_interval_milp_verified as T
rng = np.random.default_rng(11)
for K, n in ((5, 40), (6, 20), (7, 6)):
    E = list(itertools.combinations(range(K), 2)); Z = np.array([[(m >> i) & 1 for i in range(K)] for m in range(1 << K)])
    T0 = time.time(); etas = []; fb = 0; rel = []
    while len(etas) < n:
        k = rng.integers(K, 3 * K); S = rng.choice(1 << K, k, replace=False); pi = np.zeros(1 << K); pi[S] = rng.dirichlet(np.ones(k))
        qq = rng.uniform(0.1, 0.4, K); pi = 0.5 * pi + 0.5 * np.prod(np.where(Z == 1, qq, 1 - qq), axis=1)
        q = Z.T @ pi
        if np.any(q >= 0.5): continue
        e = np.array([pi @ (Z[:, i] * Z[:, j]) for i, j in E]); eps = F(int(rng.integers(0, 30)), 1000)
        qr = [F(round(x, 4)).limit_denominator(10**4) for x in q]; er = [F(round(x, 4)).limit_denominator(10**4) for x in e]
        ql = [max(F(0), x - eps) for x in qr]; qu = [min(F(1, 2), x + eps) for x in qr]
        b = [min(qu[i], qu[j], er[t] + eps) for t, (i, j) in enumerate(E)]
        g = T.certified_gap(K, ql, qu, b, time_limit=20)
        if g['lower_cert']['pi'] is None: fb += 1
        etas.append(float(g['eta'])); rel.append(float(g['eta']) / max(float(g['upper']), 1e-12))
    print(f'K={K}: {n} instances; exact lower bound reconstructed in {n-fb}/{n}; certified eta: median {np.median(etas):.4f}, '
          f'max {max(etas):.4f}, zero in {sum(x == 0 for x in etas)}; median eta/U {np.median(rel):.3f}  [{time.time()-T0:.0f}s]', flush=True)
