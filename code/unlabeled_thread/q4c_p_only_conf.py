"""Q4c: (i) confidence set with Hoeffding radius on p ONLY (unlabeled law plug-in with minimal slack, since
unlabeled data are plentiful) -- does it restore coverage with useful width?  (ii) anti family degeneracy check."""
import json, warnings
import numpy as np
from core import *
from q4_empirical import datasets, classifiers, Gf, NR, TIDX, MAJ, DICT, K, zidx, build_set, min_slack, worst_best, REPS, PAIRS
warnings.filterwarnings("ignore")
FL = flip_index(K)
out = {}
for name, X, y in datasets():
    recs = []
    for rep in range(30):
        rng = np.random.default_rng(1000 + rep)
        idx = rng.permutation(len(y)); n = len(y)
        a, b, c = int(0.4 * n), int(0.5 * n), int(0.75 * n)
        tr, lab, unl, te = idx[:a], idx[a:b], idx[b:c], idx[c:]
        F = {"lab": [], "unl": [], "te": []}
        for clf in classifiers(rep):
            clf.fit(X[tr], y[tr])
            for part, ii in (("lab", lab), ("unl", unl), ("te", te)):
                F[part].append(clf.predict(X[ii]))
        F = {k: np.array(v).T for k, v in F.items()}
        Zl = (F["lab"] == y[lab][:, None]).astype(int); Zt = (F["te"] == y[te][:, None]).astype(int)
        p_hat = Zl.mean(0)
        a_hat = np.array([np.mean(F["unl"][:, i] == F["unl"][:, j]) for i, j in PAIRS])
        phi_u = np.bincount(zidx(F["unl"]), minlength=2 ** K) / len(unl)
        m_hat = (phi_u + phi_u[FL])[REPS]
        acc = Gf @ (np.bincount(zidx(Zt), minlength=2 ** K) / len(te))
        eps_p = np.sqrt(np.log(2 * K / 0.05) / (2 * len(lab)))
        t_full = min_slack(p_hat, a_hat, m_hat, "full")
        S = build_set(p_hat, a_hat, m_hat, eps_p, 0, t_full, "full")
        lows = np.array([worst_best(Gf[k], S)[0] for k in range(NR)])
        k_rob = int(np.argmax(lows))
        r = dict(acc_rob=acc[k_rob], is_dict=k_rob in DICT)
        for tag, k in (("maj", MAJ), ("rob", k_rob)):
            lo, hi = worst_best(Gf[k], S)
            r[f"cov_{tag}"] = lo - 1e-9 <= acc[k] <= hi + 1e-9; r[f"w_{tag}"] = hi - lo; r[f"lo_{tag}"] = lo; r[f"acc_{tag}"] = acc[k]
        recs.append(r)
    d = {k: float(np.mean([r[k] for r in recs])) for k in recs[0]}
    d["eps_p"] = float(eps_p)
    out[name] = d
    print(name, {k: round(v, 4) for k, v in d.items()})
# anti-family degeneracy
from q2b_regret_and_anticorr import anti_instance  # noqa
rng = np.random.default_rng(7)
deg = {}
for K_ in (3, 5):
    Z = patterns(K_); maj = (Z.sum(1) > K_ / 2).astype(float)
    w = []; tz = []
    for _ in range(200):
        pi = anti_instance(K_, rng); p = marginals(pi, K_)
        if np.any(p <= 0.5):
            continue
        r = pairwise_joint(pi, K_)
        A, bb = set_pairwise(p, r, K_)
        lo, hi = lp_range(maj, A, bb); w.append(hi - lo)
        tz.append(max(1 - p[i] - p[j] + r[(i, j)] for (i, j) in r))
    deg[f"K{K_}"] = dict(n=len(w), frac_zero_width_majority=float(np.mean(np.array(w) < 1e-9)), mean_width=float(np.mean(w)),
                         frac_all_tandem_zero=float(np.mean(np.array(tz) < 1e-9)))
print("anti family degeneracy:", deg)
out["anti_degeneracy"] = deg
json.dump(out, open("results/q4c_summary.json", "w"), indent=2)
