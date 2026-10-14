"""Q4: small empirical sanity check with sklearn bundled / synthetic binary tasks.

Per random split: train 40% | labeled val 10% | unlabeled 25% | test 25%.
5 classifiers: LogReg, GaussianNB, DecisionTree(depth 3), kNN(7), RandomForest(50 trees, depth 5).
p_hat from labeled val (Laplace-smoothed for log-odds weights; raw for constraints);
pairwise agreement a_hat and the full even-part law m_hat from the UNLABELED split.
Rule set: all 81 self-dual positive-weight threshold rules on 5 voters.
Uncertainty sets for pi = law(Z) (Z_i = 1{clf i correct}):
  PW-plug : |marg - p_hat|_inf <= t*, |agree - a_hat|_inf <= t*  (t* = smallest slack making it feasible)
  FULL-plug: |marg - p_hat|_inf <= t*, ||even(pi) - m_hat||_1 <= t*_m  (smallest feasible)
  PW-conf : Hoeffding radii (delta = 0.05, union bound)
  FULL-conf: Hoeffding radius on p, Weissman L1 radius on m_hat (16 cells)
Report: test accuracy of max-min rules vs majority / log-odds / best single (by val acc),
interval coverage of the test accuracy, interval widths.
"""
import json, time, warnings
import numpy as np
from scipy.optimize import linprog
from sklearn.datasets import load_breast_cancer, load_digits, make_classification, make_moons
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from core import patterns, flip_index, reps, enumerate_selfdual_ltfs, ltf_table

warnings.filterwarnings("ignore")
K = 5
Z = patterns(K)
FL = flip_index(K)
REPS = reps(K)
G, _ = enumerate_selfdual_ltfs(K, n_random=50000)
Gf = G.astype(float)
NR = len(G)
TIDX = {tuple(G[k]): k for k in range(NR)}
MAJ = TIDX[tuple((Z.sum(1) >= 3).astype(int))]
DICT = [TIDX[tuple(Z[:, i].astype(int))] for i in range(K)]
PAIRS = [(i, j) for i in range(K) for j in range(i + 1, K)]
AGREE = np.array([(Z[:, i] == Z[:, j]).astype(float) for i, j in PAIRS])  # pairwise agreement rows


def datasets():
    X, y = load_breast_cancer(return_X_y=True)
    yield "breast_cancer", X, y
    X, y = load_digits(return_X_y=True)
    yield "digits_even_odd", X, (y % 2 == 0).astype(int)
    X, y = make_classification(n_samples=6000, n_features=20, n_informative=4, n_redundant=8, n_repeated=0,
                               n_clusters_per_class=3, flip_y=0.05, class_sep=0.7, random_state=0)
    yield "synthetic_corr", X, y
    X, y = make_moons(n_samples=4000, noise=0.35, random_state=0)
    rng = np.random.default_rng(0)
    X = np.hstack([X, X @ rng.normal(size=(2, 4)) + 0.5 * rng.normal(size=(4000, 4))])
    yield "moons_noisy", X, y


def classifiers(seed):
    return [make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)),
            GaussianNB(),
            DecisionTreeClassifier(max_depth=3, random_state=seed),
            make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=7)),
            RandomForestClassifier(n_estimators=50, max_depth=5, random_state=seed)]


def zidx(Zmat):
    return (Zmat.astype(int) * (2 ** np.arange(K))).sum(1)


def build_set(p_hat, a_hat, m_hat, eps_p, eps_a, eps_m, kind):
    """Return (A_eq, b_eq, A_ub, b_ub, nvar) for pi (first 32 vars) [+ aux for L1]."""
    n = 2 ** K
    if kind == "pw":
        nv = n
        A_eq = np.ones((1, nv)); b_eq = [1.0]
        A_ub = np.vstack([Z.T, -Z.T, AGREE, -AGREE])
        b_ub = np.concatenate([p_hat + eps_p, -(p_hat - eps_p), a_hat + eps_a, -(a_hat - eps_a)])
        return A_eq, np.array(b_eq), A_ub, b_ub, nv
    # full: aux u_v >= |pi(v)+pi(vbar) - m(v)| for v in REPS, sum u <= eps_m
    nr = len(REPS)
    nv = n + nr
    E = np.zeros((nr, n))
    for k, v in enumerate(REPS):
        E[k, v] = 1; E[k, FL[v]] = 1
    A_eq = np.concatenate([np.ones(n), np.zeros(nr)])[None, :]
    Zt = np.hstack([Z.T, np.zeros((K, nr))])
    A_ub = np.vstack([Zt, -Zt,
                      np.hstack([E, -np.eye(nr)]), np.hstack([-E, -np.eye(nr)]),
                      np.concatenate([np.zeros(n), np.ones(nr)])[None, :]])
    b_ub = np.concatenate([p_hat + eps_p, -(p_hat - eps_p), m_hat, -m_hat, [eps_m]])
    return A_eq, np.array([1.0]), A_ub, b_ub, nv


def min_slack(p_hat, a_hat, m_hat, kind):
    """Smallest common slack t making the set feasible (eps_p = eps_a = t, or eps_p = t, eps_m = t)."""
    lo, hi = 0.0, 1.0
    for _ in range(30):
        mid = (lo + hi) / 2
        A_eq, b_eq, A_ub, b_ub, nv = build_set(p_hat, a_hat, m_hat, mid, mid, mid, kind)
        res = linprog(np.zeros(nv), A_eq=A_eq, b_eq=b_eq, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
        if res.status == 0:
            hi = mid
        else:
            lo = mid
    return hi * 1.02 + 1e-6   # small safety margin so HiGHS reliably finds the set non-empty


def worst_best(c, S):
    A_eq, b_eq, A_ub, b_ub, nv = S
    cc = np.zeros(nv); cc[:len(c)] = c
    r1 = linprog(cc, A_eq=A_eq, b_eq=b_eq, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
    r2 = linprog(-cc, A_eq=A_eq, b_eq=b_eq, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
    if r1.status != 0 or r2.status != 0:
        raise RuntimeError(f"LP failed: {r1.status} {r2.status}")
    return r1.fun, -r2.fun


def run(n_rep=30, delta=0.05):
    allres = {}
    for name, X, y in datasets():
        t0 = time.time()
        recs = []
        for rep in range(n_rep):
            rng = np.random.default_rng(1000 + rep)
            idx = rng.permutation(len(y))
            n = len(y)
            a, b, c = int(0.4 * n), int(0.5 * n), int(0.75 * n)
            tr, lab, unl, te = idx[:a], idx[a:b], idx[b:c], idx[c:]
            clfs = classifiers(rep)
            P = {}
            for part, ii in (("lab", lab), ("unl", unl), ("te", te)):
                P[part] = []
            for clf in clfs:
                clf.fit(X[tr], y[tr])
                for part, ii in (("lab", lab), ("unl", unl), ("te", te)):
                    P[part].append(clf.predict(X[ii]))
            F = {k: np.array(v).T for k, v in P.items()}  # n x K predictions
            Zl = (F["lab"] == y[lab][:, None]).astype(int)
            Zt = (F["te"] == y[te][:, None]).astype(int)
            n_lab, n_unl = len(lab), len(unl)
            p_hat = Zl.mean(0)
            p_s = (Zl.sum(0) + 1) / (n_lab + 2)
            Fu = F["unl"]
            a_hat = np.array([np.mean(Fu[:, i] == Fu[:, j]) for i, j in PAIRS])
            phi_u = np.bincount(zidx(Fu), minlength=2 ** K) / n_unl
            m_hat = (phi_u + phi_u[FL])[REPS]
            # test-set truth
            pi_te = np.bincount(zidx(Zt), minlength=2 ** K) / len(te)
            acc_te = Gf @ pi_te
            # sets
            eps_p = np.sqrt(np.log(2 * K / (delta / 2)) / (2 * n_lab))
            eps_a = np.sqrt(np.log(2 * len(PAIRS) / (delta / 2)) / (2 * n_unl))
            eps_m = np.sqrt(2 * (np.log(2 ** 16 - 2) + np.log(2 / (delta / 2))) / n_unl)
            sets = {}
            t_pw = min_slack(p_hat, a_hat, m_hat, "pw")
            sets["pw_plug"] = build_set(p_hat, a_hat, m_hat, t_pw, t_pw, 0, "pw")
            t_full = min_slack(p_hat, a_hat, m_hat, "full")
            sets["full_plug"] = build_set(p_hat, a_hat, m_hat, t_full, 0, t_full, "full")
            sets["pw_conf"] = build_set(p_hat, a_hat, m_hat, eps_p, eps_a, 0, "pw")
            sets["full_conf"] = build_set(p_hat, a_hat, m_hat, eps_p, 0, eps_m, "full")
            rec = dict(rep=rep, p_hat=p_hat.tolist(), p_test=Zt.mean(0).tolist(), t_pw=t_pw, t_full=t_full,
                       eps_p=eps_p, eps_a=eps_a, eps_m=eps_m, acc_te=acc_te.tolist())
            k_lo = TIDX[tuple(ltf_table(np.log(p_s / (1 - p_s)), K))]
            k_best = DICT[int(np.argmax(p_hat + 1e-9 * rng.random(K)))]
            rec.update(k_lo=k_lo, k_best=k_best, k_maj=MAJ, k_oracle=int(np.argmax(acc_te)))
            for sname, S in sets.items():
                lows = np.array([worst_best(Gf[k], S)[0] for k in range(NR)])
                k_rob = int(np.argmax(lows))
                rec[f"k_rob_{sname}"] = k_rob
                rec[f"lo_rob_{sname}"] = float(lows[k_rob])
                for tag, k in (("maj", MAJ), ("lo", k_lo), ("rob", k_rob)):
                    lo_, hi_ = worst_best(Gf[k], S)
                    rec[f"int_{tag}_{sname}"] = (float(lo_), float(hi_))
            recs.append(rec)
        allres[name] = recs
        print(f"{name}: {n_rep} reps in {time.time()-t0:.0f}s (n={len(y)}, n_lab={n_lab}, n_unl={n_unl}, n_test={len(te)})")
    return allres


def summarize(allres):
    S = {}
    for name, recs in allres.items():
        acc = lambda key: np.array([r["acc_te"][r[key]] for r in recs])
        d = {}
        for lab_, key in (("best_single(val)", "k_best"), ("majority", "k_maj"), ("logodds", "k_lo"),
                          ("maxmin_pw_plug", "k_rob_pw_plug"), ("maxmin_full_plug", "k_rob_full_plug"),
                          ("maxmin_pw_conf", "k_rob_pw_conf"), ("maxmin_full_conf", "k_rob_full_conf"),
                          ("oracle_rule(test)", "k_oracle")):
            v = acc(key)
            d[lab_] = dict(mean=float(v.mean()), sd=float(v.std(ddof=1)))
        d["best_individual_test_mean"] = float(np.mean([max(r["p_test"]) for r in recs]))
        d["p_test_mean"] = np.mean([r["p_test"] for r in recs], axis=0).round(4).tolist()
        cov = {}
        for sname in ("pw_plug", "full_plug", "pw_conf", "full_conf"):
            for tag, key in (("maj", "k_maj"), ("lo", "k_lo"), ("rob", f"k_rob_{sname}")):
                ins = [r[f"int_{tag}_{sname}"][0] - 1e-9 <= r["acc_te"][r[key]] <= r[f"int_{tag}_{sname}"][1] + 1e-9 for r in recs]
                wid = [r[f"int_{tag}_{sname}"][1] - r[f"int_{tag}_{sname}"][0] for r in recs]
                lowgap = [r["acc_te"][r[key]] - r[f"int_{tag}_{sname}"][0] for r in recs]
                cov[f"{tag}_{sname}"] = dict(coverage=float(np.mean(ins)), mean_width=float(np.mean(wid)),
                                             mean_test_minus_lower=float(np.mean(lowgap)))
        d["intervals"] = cov
        d["t_pw_mean"] = float(np.mean([r["t_pw"] for r in recs]))
        d["t_full_mean"] = float(np.mean([r["t_full"] for r in recs]))
        d["eps"] = dict(eps_p=recs[0]["eps_p"], eps_a=recs[0]["eps_a"], eps_m=recs[0]["eps_m"])
        d["frac_maxmin_pw_plug_is_dictator"] = float(np.mean([r["k_rob_pw_plug"] in DICT for r in recs]))
        d["frac_maxmin_full_plug_is_dictator"] = float(np.mean([r["k_rob_full_plug"] in DICT for r in recs]))
        d["frac_logodds_is_majority"] = float(np.mean([r["k_lo"] == MAJ for r in recs]))
        S[name] = d
    return S


if __name__ == "__main__":
    allres = run()
    json.dump(allres, open("results/q4_records.json", "w"))
    S = summarize(allres)
    json.dump(S, open("results/q4_summary.json", "w"), indent=2)
    for name, d in S.items():
        print(f"\n=== {name} === p_test mean {d['p_test_mean']}; slack t_pw {d['t_pw_mean']:.4f} t_full {d['t_full_mean']:.4f}; eps {d['eps']}")
        for k in ("best_single(val)", "majority", "logodds", "maxmin_pw_plug", "maxmin_full_plug", "maxmin_pw_conf", "maxmin_full_conf", "oracle_rule(test)"):
            print(f"   {k:20s} {d[k]['mean']:.4f} +- {d[k]['sd']:.4f}")
        print(f"   dictator freq: pw_plug {d['frac_maxmin_pw_plug_is_dictator']:.2f}, full_plug {d['frac_maxmin_full_plug_is_dictator']:.2f}; logodds==majority {d['frac_logodds_is_majority']:.2f}")
        for k, v in d["intervals"].items():
            print(f"   interval {k:16s} coverage {v['coverage']:.2f}  width {v['mean_width']:.4f}  test-lower {v['mean_test_minus_lower']:.4f}")
