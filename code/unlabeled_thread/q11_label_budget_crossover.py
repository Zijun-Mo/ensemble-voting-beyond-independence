"""Q11 (completion session): how many marginal labels per classifier are needed before the max-min rule on the
confidence set stops being a dictator?  Real data (PMLB), K=5, split rep 0 as in q8; TEST set treated as the
population (exact p and m), unlabeled law treated as known (eps_m = 0), Hoeffding radius eps_p on p.
For eps_p on a grid: max-min rule over the 81 self-dual LTFs, its certified value, and its true (test) accuracy.
n_lab(eps) = log(2K/(delta/2)) / (2 eps^2)  with delta = 0.05.
"""
import json, warnings
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from core import patterns, flip_index, reps, enumerate_selfdual_ltfs, canonical_weights

warnings.filterwarnings("ignore")
import os; DATA = os.environ.get("PMLB_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_real"))
K = 5; Z = patterns(K); FL = flip_index(K); R = reps(K); N = 2 ** K
G, _ = enumerate_selfdual_ltfs(K, n_random=50000); Gf = G.astype(float)
TYPES = [str(tuple(sorted(canonical_weights(G[k], K), reverse=True))) for k in range(len(G))]
nr = len(R); P = np.zeros((nr, N))
for k, v in enumerate(R):
    P[k, v] += 1; P[k, FL[v]] += 1


def models(seed):
    return [make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)), GaussianNB(),
            DecisionTreeClassifier(max_depth=3, random_state=seed),
            make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=7)),
            RandomForestClassifier(n_estimators=50, max_depth=5, random_state=seed, n_jobs=1)]


def lowers(p, m, eps):
    A_ub = np.vstack([Z.T, -Z.T]); b_ub = np.concatenate([p + eps, -(p - eps)])
    A_eq = P; b_eq = m
    return np.array([linprog(g, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs").fun for g in Gf])


out = {}
eps_grid = [0.0, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.137]
for name in ["adult", "magic", "spambase", "phoneme", "twonorm", "churn"]:
    d = pd.read_csv(f"{DATA}/{name}.tsv.gz", sep="\t")
    y = d["target"].to_numpy(); X = d.drop(columns=["target"]).to_numpy(float); y = (y == np.unique(y)[-1]).astype(int)
    rng = np.random.default_rng(7)
    idx = rng.permutation(len(y)); n_tr = min(4000, len(y) // 2)
    tr, te = idx[:n_tr], idx[n_tr:]
    blocks = np.array_split(rng.permutation(X.shape[1]), K)
    for ens in ("same", "views"):
        F = []
        for i, clf in enumerate(models(0)):
            cols = np.arange(X.shape[1]) if ens == "same" else blocks[i]
            clf.fit(X[tr][:, cols], y[tr]); F.append(clf.predict(X[te][:, cols]))
        F = np.array(F).T; Zt = F == y[te][:, None]
        pi = np.bincount((Zt.astype(int) * 2 ** np.arange(K)).sum(1), minlength=N) / len(te)
        p = Z.T @ pi; m = (pi + pi[FL])[R]; acc = Gf @ pi
        rows = []
        for eps in eps_grid:
            L = lowers(p, m, eps); k = int(np.argmax(L))
            rows.append(dict(eps=eps, n_lab=(np.log(2 * K / 0.025) / (2 * eps ** 2)) if eps > 0 else None,
                             type=TYPES[k], cert=float(L[k]), acc=float(acc[k])))
        out[f"{name}/{ens}"] = dict(best_single=float(p.max()), oracle=float(acc.max()), oracle_type=TYPES[int(np.argmax(acc))], grid=rows)
        nd = [r for r in rows if r["type"] != "(1, 0, 0, 0, 0)"]
        print(name, ens, "best1=%.4f oracle=%.4f" % (p.max(), acc.max()),
              " | ".join(f"eps={r['eps']}: {r['type']} acc={r['acc']:.4f} cert={r['cert']:.4f}" for r in rows if r['eps'] in (0.0, 0.01, 0.03, 0.06)), flush=True)
json.dump(out, open("results/q11_results.json", "w"), indent=1)
