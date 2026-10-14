"""Q4b: separate partial-identification width from sampling error on the same real-data splits as q4.

Uses the SAME splits / classifiers as q4_empirical.py.  On the TEST split (treated as the population):
  - exact p (test labels) and exact m (test predictions) -> sharp identification intervals for
    majority / log-odds / best-single rules (pairwise info and full unlabeled info);
  - max-min rule chosen with exact (p, r) and with exact (p, m): its test accuracy.
Also: plug-in intervals widened by a 'sampling' radius is what q4 'conf' sets do; here we show the
identification-only width for reference.
"""
import json, warnings
import numpy as np
from core import patterns, flip_index, reps, set_pairwise, set_even, lp_range, lp_min, ltf_table, pairwise_joint
from q4_empirical import datasets, classifiers, Gf, NR, TIDX, MAJ, DICT, K, zidx

warnings.filterwarnings("ignore")
FL = flip_index(K)
out = {}
for name, X, y in datasets():
    rows = []
    for rep in range(30):
        rng = np.random.default_rng(1000 + rep)
        idx = rng.permutation(len(y)); n = len(y)
        a, b, c = int(0.4 * n), int(0.5 * n), int(0.75 * n)
        tr, te = idx[:a], idx[c:]
        F = []
        for clf in classifiers(rep):
            clf.fit(X[tr], y[tr]); F.append(clf.predict(X[te]))
        F = np.array(F).T
        Zt = (F == y[te][:, None]).astype(int)
        pi = np.bincount(zidx(Zt), minlength=2 ** K) / len(te)
        p = Zt.mean(0); r = pairwise_joint(pi, K); m = pi + pi[FL]
        A_p, b_p = set_pairwise(p, r, K); A_e, b_e = set_even(p, m, K)
        acc = Gf @ pi
        pc = np.clip(p, 1e-3, 1 - 1e-3)
        k_lo = TIDX[tuple(ltf_table(np.log(pc / (1 - pc)), K))]
        k_best = DICT[int(np.argmax(p))]
        lows_p = np.array([lp_min(Gf[k], A_p, b_p) for k in range(NR)])
        lows_e = np.array([lp_min(Gf[k], A_e, b_e) for k in range(NR)])
        rec = dict(acc_maj=acc[MAJ], acc_lo=acc[k_lo], acc_best=acc[k_best], acc_oracle=acc.max(),
                   acc_mm_pw=acc[int(np.argmax(lows_p))], acc_mm_full=acc[int(np.argmax(lows_e))])
        for tag, k in (("maj", MAJ), ("lo", k_lo)):
            lo1, hi1 = lp_range(Gf[k], A_p, b_p); lo2, hi2 = lp_range(Gf[k], A_e, b_e)
            rec[f"w_{tag}_pw"] = hi1 - lo1; rec[f"w_{tag}_full"] = hi2 - lo2
            rec[f"lowgap_{tag}_pw"] = acc[k] - lo1; rec[f"lowgap_{tag}_full"] = acc[k] - lo2
        rows.append(rec)
    d = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
    out[name] = d
    print(name, {k: round(v, 4) for k, v in d.items()})
json.dump(out, open("results/q4b_summary.json", "w"), indent=2)
