"""Q8b (completion session): which labeled-information regime makes partial identification worthwhile?
Same datasets / ensembles as q8 ('same' and 'views'), K=5, 81 self-dual weighted rules.
Regimes (labeled budget n_lab per classifier):
  MARG  : marginal-only labels (each classifier's accuracy measured on its OWN labeled sample of size n_lab,
          e.g. vendor/leaderboard numbers) + unlabeled law -> confidence set C(eps_p, eps_m); max-min rule.
  JOINT : one shared labeled sample of size n_lab (joint correctness patterns observed) WITHOUT unlabeled data:
          per-rule Hoeffding + union bound over 81 rules; pick argmax of certified lower bound.
  HYBRID: JOINT constraints (|acc_hat_k - G_k.pi| <= eps_81, all k) + marginal box + unlabeled ball -> LP.
Report test accuracy of chosen rule, its certified lower bound, coverage (cert <= test acc), for n_lab in {200, 1000}.
"""
import json, time, warnings, sys
import numpy as np
import pandas as pd
from collections import Counter
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
G, _ = enumerate_selfdual_ltfs(K, n_random=50000); Gf = G.astype(float); NR = len(G)
TYPES = [str(tuple(sorted(canonical_weights(G[k], K), reverse=True))) for k in range(NR)]
NAMES = sys.argv[1].split(",")
NREP = int(sys.argv[2])
DELTA = 0.05
nr = len(R); P = np.zeros((nr, N))
for k, v in enumerate(R):
    P[k, v] += 1; P[k, FL[v]] += 1


def models(seed):
    return [make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)), GaussianNB(),
            DecisionTreeClassifier(max_depth=3, random_state=seed),
            make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=7)),
            RandomForestClassifier(n_estimators=50, max_depth=5, random_state=seed, n_jobs=1)]


def zidx(M):
    return (M.astype(int) * (2 ** np.arange(K))).sum(1)


def lower_all(p_hat, eps_p, m_hat, eps_m, acc_hat=None, eps_r=None):
    rows = [np.hstack([Z.T, np.zeros((K, nr))]), np.hstack([-Z.T, np.zeros((K, nr))])]
    rhs = [p_hat + eps_p, -(p_hat - eps_p)]
    if m_hat is not None:
        rows += [np.hstack([P, -np.eye(nr)]), np.hstack([-P, -np.eye(nr)]), np.concatenate([np.zeros(N), np.ones(nr)])[None, :]]
        rhs += [m_hat, -m_hat, [eps_m]]
    if acc_hat is not None:
        rows += [np.hstack([Gf, np.zeros((NR, nr))]), np.hstack([-Gf, np.zeros((NR, nr))])]
        rhs += [acc_hat + eps_r, -(acc_hat - eps_r)]
    A_ub = np.vstack(rows); b_ub = np.concatenate(rhs)
    A_eq = np.concatenate([np.ones(N), np.zeros(nr)])[None, :]
    out = []
    for g in Gf:
        res = linprog(np.concatenate([g, np.zeros(nr)]), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=(0, None), method="highs")
        out.append(res.fun if res.status == 0 else np.nan)
    return np.array(out)


def load(name):
    d = pd.read_csv(f"{DATA}/{name}.tsv.gz", sep="\t")
    y = d["target"].to_numpy(); X = d.drop(columns=["target"]).to_numpy(float)
    return X, (y == np.unique(y)[-1]).astype(int)


results = {}
for name in NAMES:
    X, y = load(name); t0 = time.time(); recs = []
    for rep in range(NREP):
        rng = np.random.default_rng(10 * rep + 7)
        idx = rng.permutation(len(y)); n_pool = 5 * 1000
        n_unl = min(5000, int(0.3 * len(y)))
        rest = len(y) - n_unl - n_pool
        n_tr = min(4000, max(rest // 2, 500))
        unl = idx[:n_unl]; tr = idx[n_unl:n_unl + n_tr]; pool = idx[n_unl + n_tr:n_unl + n_tr + n_pool]
        te = idx[n_unl + n_tr + n_pool:]
        if len(te) < 500:  # small dataset: shrink pool
            pool = idx[n_unl + n_tr:n_unl + n_tr + (len(y) - n_unl - n_tr) // 2]; te = idx[n_unl + n_tr + len(pool):]
        d = X.shape[1]; blocks = np.array_split(rng.permutation(d), K)
        for ens in ("same", "views"):
            Fp, Fu, Ft = [], [], []
            for i, clf in enumerate(models(rep)):
                cols = np.arange(d) if ens == "same" else blocks[i]
                clf.fit(X[tr][:, cols], y[tr])
                Fp.append(clf.predict(X[pool][:, cols])); Fu.append(clf.predict(X[unl][:, cols])); Ft.append(clf.predict(X[te][:, cols]))
            Fp, Fu, Ft = np.array(Fp).T, np.array(Fu).T, np.array(Ft).T
            Zp = Fp == y[pool][:, None]; Zt = Ft == y[te][:, None]
            acc = Gf @ (np.bincount(zidx(Zt), minlength=N) / len(te))
            phi = np.bincount(zidx(Fu), minlength=N) / n_unl; m_hat = (phi + phi[FL])[R]
            eps_m = np.sqrt(2 * (nr * np.log(2) + np.log(2 / DELTA)) / n_unl)
            for n_lab in (200, 1000):
                if n_lab * K > len(pool):
                    continue
                # MARG: classifier i uses its own disjoint block of n_lab labeled points
                p_marg = np.array([Zp[i * n_lab:(i + 1) * n_lab, i].mean() for i in range(K)])
                eps_p = np.sqrt(np.log(2 * K / (DELTA / 2)) / (2 * n_lab))
                L_marg = lower_all(p_marg, eps_p, m_hat, eps_m)
                # JOINT: shared sample of n_lab points
                Zj = Zp[:n_lab]; pij = np.bincount(zidx(Zj), minlength=N) / n_lab
                acc_hat = Gf @ pij
                eps_r = np.sqrt(np.log(2 * NR / DELTA) / (2 * n_lab))
                L_joint = acc_hat - eps_r
                eps_pj = np.sqrt(np.log(2 * K / (DELTA / 2)) / (2 * n_lab))
                L_hyb = lower_all(Zj.mean(0), eps_pj, m_hat, eps_m, acc_hat, np.sqrt(np.log(2 * NR / (DELTA / 2)) / (2 * n_lab)))
                rec = dict(rep=rep, ens=ens, n_lab=n_lab, oracle=acc.max(), best_test_single=acc[[TYPES.index("(1, 0, 0, 0, 0)")]].max())
                for tag, L in (("marg", L_marg), ("joint", L_joint), ("hyb", L_hyb)):
                    k = int(np.nanargmax(L))
                    rec[f"acc_{tag}"] = float(acc[k]); rec[f"cert_{tag}"] = float(L[k]); rec[f"type_{tag}"] = TYPES[k]
                    rec[f"cov_{tag}"] = bool(np.all(L[~np.isnan(L)] <= acc[~np.isnan(L)] + 1e-9))
                recs.append(rec)
    results[name] = recs
    print(f"{name} {time.time() - t0:.0f}s", flush=True)

summ = {}
for name, recs in results.items():
    for ens in ("same", "views"):
        for n_lab in (200, 1000):
            rr = [r for r in recs if r["ens"] == ens and r["n_lab"] == n_lab]
            if not rr:
                continue
            s = {k: float(np.mean([r[k] for r in rr])) for k in rr[0] if k.startswith(("acc", "cert", "cov", "oracle"))}
            for tag in ("marg", "joint", "hyb"):
                s[f"types_{tag}"] = dict(Counter(r[f"type_{tag}"] for r in rr))
            s["n"] = len(rr)
            summ[f"{name}/{ens}/n{n_lab}"] = s
            print(name, ens, n_lab, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()}, flush=True)
json.dump(summ, open(f"results/q8b_summary_{'_'.join(NAMES)}.json", "w"), indent=1)
