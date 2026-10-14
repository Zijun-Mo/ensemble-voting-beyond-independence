"""Q2c: symmetric case (p_i = p, rho_ij = rho).  Which self-dual LTF type is max-min optimal over P(p, r)?
In the symmetric case all rules of the same weight type have the same worst case (permutation invariance), so
one representative per type suffices.  K=3: analytic threshold rho*(p) = ((5p-2)/3 - p^2)/(p(1-p)).
Log-odds weights are equal here, so log-odds == simple majority everywhere."""
import json
import numpy as np
from core import patterns, set_pairwise, lp_min, ltf_table

out = {}
for K in (3, 5):
    Z = patterns(K)
    types = [(1, 0, 0), (1, 1, 1)] if K == 3 else [(1, 0, 0, 0, 0), (1, 1, 1, 0, 0), (2, 1, 1, 1, 0), (1, 1, 1, 1, 1),
                                                    (3, 1, 1, 1, 1), (2, 2, 1, 1, 1), (3, 2, 2, 1, 1)]
    tables = [ltf_table(np.array(t, float), K).astype(float) for t in types]
    ps = np.linspace(0.505, 0.995, 99)
    rhos = np.linspace(-0.5, 1.0, 151)
    best = np.full((len(rhos), len(ps)), -1)
    worst_vals = np.full((len(rhos), len(ps), len(types)), np.nan)
    for a, p in enumerate(ps):
        for b, rho in enumerate(rhos):
            r = p * p + rho * p * (1 - p)
            rr = {(i, j): r for i in range(K) for j in range(i + 1, K)}
            A, bb = set_pairwise(np.full(K, p), rr, K)
            vals = []
            ok = True
            for t in tables:
                v = lp_min(t, A, bb)
                if v is None:
                    ok = False; break
                vals.append(v)
            if not ok:
                continue
            vals = np.array(vals)
            worst_vals[b, a] = vals
            best[b, a] = int(np.argmax(vals + 1e-9 * np.arange(len(vals))[::-1] * 0))
    np.save(f"results/q2c_best_K{K}.npy", best)
    np.save(f"results/q2c_worst_K{K}.npy", worst_vals)
    feas = best >= 0
    frac = {str(types[k]): float(np.mean(best[feas] == k)) for k in range(len(types))}
    out[f"K{K}"] = dict(types=[str(t) for t in types], frac_of_feasible_grid=frac, ps=ps.tolist(), rhos=rhos.tolist())
    print(f"K={K}: fraction of feasible (p,rho) grid where each type is max-min optimal: {frac}")
    if K == 3:
        thr = ((5 * ps - 2) / 3 - ps ** 2) / (ps * (1 - ps))
        # check analytic threshold: majority (index 1) optimal iff rho < thr
        agree = []
        for a in range(len(ps)):
            for b in range(len(rhos)):
                if best[b, a] < 0 or abs(rhos[b] - thr[a]) < 1e-6:
                    continue
                agree.append((best[b, a] == 1) == (rhos[b] < thr[a]))
        out["K3_threshold_agreement"] = float(np.mean(agree))
        print(f"K=3 analytic threshold agreement on grid: {np.mean(agree):.4f}")
    if K == 5:
        # where does 5-majority (index 3) become optimal? smallest rho for p values
        for p_show in (0.6, 0.7, 0.8, 0.9):
            a = int(np.argmin(np.abs(ps - p_show)))
            col = best[:, a]
            reg = {str(types[k]): [float(rhos[i]) for i in np.where(col == k)[0][[0, -1]]] if np.any(col == k) else None
                   for k in range(len(types))}
            out[f"K5_p{p_show}_rho_ranges"] = reg
            print(f"  K=5 p={ps[a]:.3f}: rho ranges per optimal type: {reg}")
json.dump(out, open("results/q2c_results.json", "w"), indent=2)
