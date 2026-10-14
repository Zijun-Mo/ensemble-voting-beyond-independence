"""Q8 (completion session): real-data experiment (PMLB copies of UCI datasets: adult, magic, spambase,
phoneme, twonorm, churn; downloaded to data_real/ or $PMLB_DATA_DIR, read with pandas only).

Two K=5 ensembles per dataset:
  'same'  : LR, GaussianNB, Tree(d3), kNN(7), RF(50,d5) on ALL features (positively dependent errors);
  'views' : the same five model classes, each trained on its OWN random disjoint feature block (multi-view,
            weaker but more complementary voters).
Budget: n_lab = 200 labeled (estimate p), n_unl = 5000 unlabeled (estimate m), rest -> train (<= 4000) / test.
Methods (rule chosen from all 81 self-dual weighted threshold rules on 5 voters):
  best single (val), majority, log-odds (Laplace p), max-min on plug-in E(m_hat, p_hat) [min slack],
  max-min on confidence set C(eps_p, eps_m) (Hoeffding on p, Weissman L1 on m), max-min-regret-free 'oracle'.
Bounds on majority error reported with TEST-set quantities as the population (oracle comparison):
  FO 2E[L], TND 4E[LL'], C-bound, min(C, E[L]+E[D]) (Thm 5), finite-K Thm 6 (mu,s), sharp|p only,
  sharp|p+pairwise, sharp|p+full unlabeled; plus true error.
Finite-sample: coverage of the realised test accuracy by [L_conf, U_conf] for majority and for the chosen rule.
"""
import json, time, warnings, sys
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
from core import patterns, flip_index, reps, enumerate_selfdual_ltfs, ltf_table, canonical_weights, set_pairwise, set_even, lp_range

warnings.filterwarnings("ignore")
import os; DATA = os.environ.get("PMLB_DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_real"))
K = 5; Z = patterns(K); FL = flip_index(K); R = reps(K); N = 2 ** K
G, _ = enumerate_selfdual_ltfs(K, n_random=50000); Gf = G.astype(float); NR = len(G)
TIDX = {tuple(G[k]): k for k in range(NR)}
MAJ = TIDX[tuple((Z.sum(1) >= 3).astype(int))]
DICT = [TIDX[tuple(Z[:, i].astype(int))] for i in range(K)]
TYPES = [str(tuple(sorted(canonical_weights(G[k], K), reverse=True))) for k in range(NR)]
PAIRS = [(i, j) for i in range(K) for j in range(i + 1, K)]
NAMES = sys.argv[1].split(",") if len(sys.argv) > 1 else ["adult", "magic", "spambase", "phoneme", "twonorm", "churn"]
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 10


def models(seed):
    return [make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)), GaussianNB(),
            DecisionTreeClassifier(max_depth=3, random_state=seed),
            make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=7)),
            RandomForestClassifier(n_estimators=50, max_depth=5, random_state=seed, n_jobs=1)]


def zidx(M):
    return (M.astype(int) * (2 ** np.arange(K))).sum(1)


def conf_set(p_hat, m_hat, eps_p, eps_m):
    nr = len(R); P = np.zeros((nr, N))
    for k, v in enumerate(R):
        P[k, v] += 1; P[k, FL[v]] += 1
    A_ub = np.vstack([np.hstack([Z.T, np.zeros((K, nr))]), np.hstack([-Z.T, np.zeros((K, nr))]),
                      np.hstack([P, -np.eye(nr)]), np.hstack([-P, -np.eye(nr)]),
                      np.concatenate([np.zeros(N), np.ones(nr)])[None, :]])
    b_ub = np.concatenate([p_hat + eps_p, -(p_hat - eps_p), m_hat, -m_hat, [eps_m]])
    A_eq = np.concatenate([np.ones(N), np.zeros(nr)])[None, :]
    return A_ub, b_ub, A_eq, N + nr


def rng_on(g, S):
    A_ub, b_ub, A_eq, nv = S
    c = np.concatenate([g, np.zeros(nv - N)])
    lo = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=(0, None), method="highs")
    hi = linprog(-c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=(0, None), method="highs")
    if lo.status or hi.status:
        return np.nan, np.nan
    return lo.fun, -hi.fun


def feasible(p_hat, m_hat, t):
    A_ub, b_ub, A_eq, nv = conf_set(p_hat, m_hat, t, t)
    return linprog(np.zeros(nv), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=(0, None), method="highs").status == 0


def min_slack(p_hat, m_hat):
    lo, hi = 0.0, 1.0
    for _ in range(25):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if feasible(p_hat, m_hat, mid) else (mid, hi)
    return hi * 1.02 + 1e-6


def thm6(mu, s):
    r = (K + 1) // 2; m1, m2 = K * mu, K * K * s
    c = [1.0, ((r + K) * m1 - m2) / (r * K)]
    for j in range(r - 1):
        c.append((m2 - (2 * j + 1) * m1 + j * (j + 1)) / ((r - j) * (r - j - 1)))
    return min(c)


def oracle_bounds(Zt):
    """uniform-weight majority error bounds from test-set (population) quantities."""
    L = 1 - Zt.mean(0); mu = L.mean()
    E = 1 - Zt
    tand = (E.T @ E) / len(E)              # K x K, diagonal = L_i
    s = tand.mean()                        # E_{rho^2} L(h,h') incl. diagonal
    Dbar = 2 * (mu - s)                    # E_{rho^2}[D]
    C = (s - mu ** 2) / (s - mu + 0.25) if s - mu + 0.25 > 0 else 1.0
    pi = np.bincount(zidx(Zt), minlength=N) / len(Zt); p = 1 - L
    maj = (Z.sum(1) >= 3).astype(float)
    A, b = np.vstack([np.ones(N), Z.T]), np.concatenate([[1.0], p])
    lo_p = linprog(maj, A_eq=A, b_eq=b, bounds=(0, None), method="highs").fun
    r = {(i, j): float(pi @ (Z[:, i] * Z[:, j])) for i, j in PAIRS}
    A, b = set_pairwise(p, r, K); lo_pw = linprog(maj, A_eq=A, b_eq=b, bounds=(0, None), method="highs").fun
    A, b = set_even(p, pi + pi[FL], K); lo_full = linprog(maj, A_eq=A, b_eq=b, bounds=(0, None), method="highs").fun
    return dict(true=float(1 - maj @ pi), FO=min(1, 2 * mu), TND=min(1, 4 * s), C=min(1, C), CminLin=min(1, C, mu + Dbar),
                thm6=thm6(mu, s), sharp_p=1 - lo_p, sharp_pw=1 - lo_pw, sharp_full=1 - lo_full, gibbs=mu, dis=Dbar)


def load(name):
    d = pd.read_csv(f"{DATA}/{name}.tsv.gz", sep="\t")
    y = d["target"].to_numpy(); X = d.drop(columns=["target"]).to_numpy(float)
    vals = np.unique(y)
    y = (y == vals[-1]).astype(int)
    return X, y


results = {}
for name in NAMES:
    X, y = load(name)
    t0 = time.time()
    recs = []
    for rep in range(NREP):
        rng = np.random.default_rng(10 * rep + 7)
        idx = rng.permutation(len(y))
        n_lab = 200; n_unl = min(5000, int(0.4 * (len(y) - n_lab)))
        n_rest = len(y) - n_lab - n_unl
        n_tr = min(4000, n_rest // 2)
        lab, unl, tr, te = idx[:n_lab], idx[n_lab:n_lab + n_unl], idx[n_lab + n_unl:n_lab + n_unl + n_tr], idx[n_lab + n_unl + n_tr:]
        d = X.shape[1]
        perm = rng.permutation(d)
        blocks = np.array_split(perm, K) if d >= K else [np.arange(d)] * K
        for ens in ("same", "views"):
            F = {"lab": [], "unl": [], "te": []}
            for i, clf in enumerate(models(rep)):
                cols = np.arange(d) if ens == "same" else blocks[i]
                clf.fit(X[tr][:, cols], y[tr])
                for part, ii in (("lab", lab), ("unl", unl), ("te", te)):
                    F[part].append(clf.predict(X[ii][:, cols]))
            F = {k: np.array(v).T for k, v in F.items()}
            Zl = (F["lab"] == y[lab][:, None]); Zt = (F["te"] == y[te][:, None])
            p_hat = Zl.mean(0); p_s = (Zl.sum(0) + 1) / (n_lab + 2)
            phi_u = np.bincount(zidx(F["unl"]), minlength=N) / n_unl
            m_hat = (phi_u + phi_u[FL])[R]
            pi_te = np.bincount(zidx(Zt), minlength=N) / len(te)
            acc = Gf @ pi_te
            eps_p = np.sqrt(np.log(2 * K / 0.025) / (2 * n_lab))
            eps_m = np.sqrt(2 * (len(R) * np.log(2) + np.log(1 / 0.025)) / n_unl)
            t_plug = min_slack(p_hat, m_hat)
            S_plug = conf_set(p_hat, m_hat, t_plug, t_plug)
            S_conf = conf_set(p_hat, m_hat, eps_p, max(eps_m, t_plug))
            lo_plug = np.array([rng_on(g, S_plug)[0] for g in Gf])
            ints_conf = np.array([rng_on(g, S_conf) for g in Gf])
            k_lo = TIDX[tuple(ltf_table(np.log(p_s / (1 - p_s)), K))]
            k_best = DICT[int(np.argmax(p_hat + 1e-9 * rng.random(K)))]
            k_mmp = int(np.nanargmax(lo_plug)); k_mmc = int(np.nanargmax(ints_conf[:, 0]))
            rec = dict(rep=rep, ens=ens, p_test=Zt.mean(0).tolist(), eps_p=eps_p, eps_m=eps_m, t_plug=t_plug,
                       acc_best=acc[k_best], acc_maj=acc[MAJ], acc_lo=acc[k_lo], acc_mmp=acc[k_mmp], acc_mmc=acc[k_mmc],
                       acc_oracle=acc.max(), type_mmp=TYPES[k_mmp], type_mmc=TYPES[k_mmc], type_lo=TYPES[k_lo],
                       type_oracle=TYPES[int(np.argmax(acc))],
                       cert_mmc=float(ints_conf[k_mmc, 0]), cert_maj=float(ints_conf[MAJ, 0]),
                       cert_best=float(ints_conf[k_best, 0]),
                       cov_maj=bool(ints_conf[MAJ, 0] - 1e-9 <= acc[MAJ] <= ints_conf[MAJ, 1] + 1e-9),
                       cov_mmc=bool(ints_conf[k_mmc, 0] - 1e-9 <= acc[k_mmc] <= ints_conf[k_mmc, 1] + 1e-9),
                       cov_all=float(np.mean((ints_conf[:, 0] - 1e-9 <= acc) & (acc <= ints_conf[:, 1] + 1e-9))),
                       w_maj=float(ints_conf[MAJ, 1] - ints_conf[MAJ, 0]))
            rec.update({f"ob_{k}": v for k, v in oracle_bounds(Zt.astype(float)).items()})
            recs.append(rec)
    results[name] = recs
    print(f"{name}: {time.time() - t0:.0f}s, n={len(y)}, d={X.shape[1]}", flush=True)
    json.dump(results, open(f"results/q8_records_{'_'.join(NAMES)}.json", "w"), default=float)

# summary
summ = {}
for name, recs in results.items():
    for ens in ("same", "views"):
        rr = [r for r in recs if r["ens"] == ens]
        keys = [k for k in rr[0] if isinstance(rr[0][k], (float, int, bool, np.floating, np.bool_)) and k != "rep"]
        s = {k: float(np.mean([r[k] for r in rr])) for k in keys}
        s.update({f"sd_{k}": float(np.std([r[k] for r in rr])) for k in ("acc_best", "acc_maj", "acc_lo", "acc_mmp", "acc_mmc")})
        from collections import Counter
        s["types_mmc"] = dict(Counter(r["type_mmc"] for r in rr)); s["types_oracle"] = dict(Counter(r["type_oracle"] for r in rr))
        summ[f"{name}/{ens}"] = s
        print(name, ens, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()
                          if k.startswith(("acc", "cert", "cov", "w_", "ob_", "types"))}, flush=True)
json.dump(summ, open(f"results/q8_summary_{'_'.join(NAMES)}.json", "w"), indent=1)
