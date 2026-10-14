"""Q2: robust (max-min) weighted votes from (p, pairwise agreement) [and from the full unlabeled law].

For K=3 and K=5, on random instances from three families, compare:
  robust_pw   : argmax over self-dual LTFs of worst-case accuracy over P(p, r)          (pairwise info)
  robust_full : argmax over self-dual LTFs of worst-case accuracy over E(m, p)          (full unlabeled law)
  logodds     : weights log(p/(1-p))  (Nitzan-Paroush; optimal under independence)
  majority    : equal weights
  best        : best single classifier (argmax p)
Reported: worst-case accuracy (pairwise set), true accuracy, how often rules differ,
plus the value of randomised max-min (mixtures of LTFs / randomised self-dual g / Balsubramani-Freund game).
K=3: verify closed form 'robust picks majority iff D > e_(2)+e_(3)' (D = P(not unanimous), e sorted errors asc).
"""
import json, time
import numpy as np
from core import *

rng = np.random.default_rng(2026)


def family_instance(K, fam, rng):
    Z = patterns(K)
    if fam == "independent":
        p = rng.uniform(0.55, 0.9, size=K)
        pi = np.prod(np.where(Z == 1, p, 1 - p), axis=1)
    elif fam == "mixture":
        pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3), acc_range=(0.55, 0.9))
    elif fam == "bloc":
        # first nb classifiers share a latent correct/incorrect draw with prob c; others independent
        nb = 2 if K == 3 else 3
        c = rng.uniform(0.3, 0.95)
        pb = rng.uniform(0.6, 0.85)
        po = rng.uniform(0.55, 0.85, size=K - nb)
        pi = np.zeros(2 ** K)
        for zb in (0, 1):  # common draw
            pc = pb if zb == 1 else 1 - pb
            q = np.concatenate([np.full(nb, c * zb + (1 - c) * pb), po])
            pi += pc * np.prod(np.where(Z == 1, q, 1 - q), axis=1)
    pi = np.clip(pi, 0, None)
    return pi / pi.sum()


def run(K, n_per_family=300):
    G, W = enumerate_selfdual_ltfs(K, n_random=50000 if K == 5 else 2000)
    Z = patterns(K)
    S = 2 * Z - 1
    nr = len(G)
    maj_idx = [k for k in range(nr) if np.array_equal(G[k], (Z.sum(1) > K / 2).astype(int))][0]
    dict_idx = [[k for k in range(nr) if np.array_equal(G[k], Z[:, i].astype(int))][0] for i in range(K)]
    table_index = {tuple(G[k]): k for k in range(nr)}
    print(f"K={K}: {nr} self-dual positive-weight threshold rules")
    recs = []
    for fam in ("independent", "mixture", "bloc"):
        for _ in range(n_per_family):
            pi = family_instance(K, fam, rng)
            p = marginals(pi, K)
            if np.any(p <= 0.5):
                continue
            r = pairwise_joint(pi, K)
            m = pi + pi[flip_index(K)]
            A_p, b_p = set_pairwise(p, r, K)
            A_e, b_e = set_even(p, m, K)
            true = G @ pi
            lo_pw = np.array([lp_min(G[k].astype(float), A_p, b_p) for k in range(nr)])
            lo_full = np.array([lp_min(G[k].astype(float), A_e, b_e) for k in range(nr)]) if K > 3 else lo_pw
            w_lo = np.log(p / (1 - p))
            k_lo = table_index[tuple(ltf_table(w_lo, K))]
            k_rpw = int(np.argmax(lo_pw + 1e-12 * true * 0))  # ties broken by first index
            k_rfull = int(np.argmax(lo_full))
            k_best = dict_idx[int(np.argmax(p))]
            mix_pw, _ = maxmin_mixture(G.astype(float), A_p, b_p)
            rsd_pw, _ = maxmin_selfdual_randomized(A_p, b_p, K)
            phi = 0.5 * m  # balanced symmetric representation; BF set depends only on (m,p) for prior unknown
            bf, _ = bf_value(phi, p, K)
            rec = dict(fam=fam, p=p.tolist(),
                       D=float(1 - m[0]) if K == 3 else None,
                       k_rpw=k_rpw, k_rfull=k_rfull, k_lo=k_lo, k_maj=maj_idx, k_best=k_best,
                       true=true.tolist(), lo_pw=lo_pw.tolist(), lo_full=lo_full.tolist(),
                       mix_pw=mix_pw, rsd_pw=rsd_pw, bf_full=bf,
                       oracle_true=float(true.max()))
            recs.append(rec)
    return G, W, recs


def summarize(K, G, recs):
    out = {}
    for fam in ("independent", "mixture", "bloc", "all"):
        R = [r for r in recs if fam == "all" or r["fam"] == fam]
        if not R:
            continue
        def col(key_k, arr):
            return np.array([r[arr][r[key_k]] for r in R])
        s = dict(n=len(R))
        for name, key in (("robust_pw", "k_rpw"), ("robust_full", "k_rfull"), ("logodds", "k_lo"),
                          ("majority", "k_maj"), ("best_single", "k_best")):
            s[f"{name}_true_mean"] = float(col(key, "true").mean())
            s[f"{name}_worst_pw_mean"] = float(col(key, "lo_pw").mean())
            s[f"{name}_worst_full_mean"] = float(col(key, "lo_full").mean())
            s[f"{name}_regret_vs_oracle_true_mean"] = float(np.mean([r["oracle_true"] - r["true"][r[key]] for r in R]))
        s["frac_robust_pw_ne_logodds"] = float(np.mean([r["k_rpw"] != r["k_lo"] for r in R]))
        s["frac_robust_full_ne_logodds"] = float(np.mean([r["k_rfull"] != r["k_lo"] for r in R]))
        s["frac_robust_pw_ne_robust_full"] = float(np.mean([r["k_rpw"] != r["k_rfull"] for r in R]))
        s["frac_robust_pw_is_dictator"] = float(np.mean([G[r["k_rpw"]].sum() == 2 ** (K - 1) and
                                                         any(np.array_equal(G[r["k_rpw"]], patterns(K)[:, i].astype(int)) for i in range(K)) for r in R]))
        s["frac_logodds_is_dictator"] = float(np.mean([any(np.array_equal(G[r["k_lo"]], patterns(K)[:, i].astype(int)) for i in range(K)) for r in R]))
        s["robust_pw_true_minus_logodds_true_mean"] = float(np.mean([r["true"][r["k_rpw"]] - r["true"][r["k_lo"]] for r in R]))
        s["robust_pw_worst_minus_logodds_worst_mean"] = float(np.mean([r["lo_pw"][r["k_rpw"]] - r["lo_pw"][r["k_lo"]] for r in R]))
        diff = [r for r in R if r["k_rpw"] != r["k_lo"]]
        if diff:
            s["when_differ_true_gain_mean"] = float(np.mean([r["true"][r["k_rpw"]] - r["true"][r["k_lo"]] for r in diff]))
            s["when_differ_frac_robust_true_better"] = float(np.mean([r["true"][r["k_rpw"]] > r["true"][r["k_lo"]] + 1e-12 for r in diff]))
            s["when_differ_worst_gain_mean"] = float(np.mean([r["lo_pw"][r["k_rpw"]] - r["lo_pw"][r["k_lo"]] for r in diff]))
        s["mixture_gain_over_det_mean"] = float(np.mean([r["mix_pw"] - max(r["lo_pw"]) for r in R]))
        s["mixture_gain_over_det_max"] = float(np.max([r["mix_pw"] - max(r["lo_pw"]) for r in R]))
        s["randomized_selfdual_gain_over_det_mean"] = float(np.mean([r["rsd_pw"] - max(r["lo_pw"]) for r in R]))
        s["randomized_selfdual_gain_over_det_max"] = float(np.max([r["rsd_pw"] - max(r["lo_pw"]) for r in R]))
        s["bf_full_gain_over_det_full_mean"] = float(np.mean([r["bf_full"] - max(r["lo_full"]) for r in R]))
        s["bf_full_gain_over_det_full_max"] = float(np.max([r["bf_full"] - max(r["lo_full"]) for r in R]))
        out[fam] = s
    return out


results = {}
for K in (3, 5):
    t0 = time.time()
    G, W, recs = run(K, n_per_family=300 if K == 3 else 200)
    summ = summarize(K, G, recs)
    results[f"K{K}"] = summ
    print(f"K={K} done in {time.time()-t0:.1f}s")
    for fam, s in summ.items():
        print(f"  [{fam}] n={s['n']}: robust_pw!=logodds {s['frac_robust_pw_ne_logodds']:.3f}; "
              f"true acc: robust_pw {s['robust_pw_true_mean']:.4f}, robust_full {s['robust_full_true_mean']:.4f}, "
              f"logodds {s['logodds_true_mean']:.4f}, maj {s['majority_true_mean']:.4f}, best {s['best_single_true_mean']:.4f}; "
              f"worst(pw): robust {s['robust_pw_worst_pw_mean']:.4f}, logodds {s['logodds_worst_pw_mean']:.4f}, "
              f"maj {s['majority_worst_pw_mean']:.4f}, best {s['best_single_worst_pw_mean']:.4f}")
        print(f"       mixture gain mean/max {s['mixture_gain_over_det_mean']:.4f}/{s['mixture_gain_over_det_max']:.4f}; "
              f"rand self-dual gain {s['randomized_selfdual_gain_over_det_mean']:.4f}/{s['randomized_selfdual_gain_over_det_max']:.4f}; "
              f"BF(full) gain {s['bf_full_gain_over_det_full_mean']:.4f}/{s['bf_full_gain_over_det_full_max']:.4f}; "
              f"robust_dict {s['frac_robust_pw_is_dictator']:.3f}, logodds_dict {s['frac_logodds_is_dictator']:.3f}")
        if "when_differ_true_gain_mean" in s:
            print(f"       when differ: true gain {s['when_differ_true_gain_mean']:.4f}, robust better in {s['when_differ_frac_robust_true_better']:.3f}, worst gain {s['when_differ_worst_gain_mean']:.4f}")
    if K == 3:
        # closed-form check
        ok = 0; n = 0
        for r in recs:
            p = np.array(r["p"]); e = np.sort(1 - p)  # ascending errors
            pred_maj = r["D"] > e[1] + e[2] + 1e-12
            lo = np.array(r["lo_pw"])
            is_maj = lo[r["k_maj"]] > max(lo[k] for k in range(len(lo)) if k != r["k_maj"]) + 1e-12
            tie = abs(r["D"] - e[1] - e[2]) < 1e-9
            if not tie:
                ok += int(pred_maj == is_maj); n += 1
        results["K3_closed_form_rule_check"] = dict(n=n, agree=ok)
        print(f"  K=3 closed-form choice rule agrees with LP in {ok}/{n} instances")
        # log-odds vs robust disagreement types
        cnt = dict(lo_maj_rob_dict=0, lo_dict_rob_maj=0)
        for r in recs:
            if r["k_lo"] == r["k_maj"] and r["k_rpw"] != r["k_maj"]:
                cnt["lo_maj_rob_dict"] += 1
            if r["k_lo"] != r["k_maj"] and r["k_rpw"] == r["k_maj"]:
                cnt["lo_dict_rob_maj"] += 1
        results["K3_disagreement_types"] = cnt
        print("  K=3 disagreement types:", cnt)
    if K == 5:
        # which rule types are chosen by robust_pw
        from collections import Counter
        typ = Counter()
        for r in recs:
            w = canonical_weights(G[r["k_rpw"]], K, maxw=5)
            typ[tuple(sorted(w, reverse=True)) if w else None] += 1
        results["K5_robust_pw_rule_types"] = {str(k): v for k, v in typ.items()}
        typ2 = Counter()
        for r in recs:
            w = canonical_weights(G[r["k_lo"]], K, maxw=5)
            typ2[tuple(sorted(w, reverse=True)) if w else None] += 1
        results["K5_logodds_rule_types"] = {str(k): v for k, v in typ2.items()}
        print("  K=5 robust_pw rule types:", dict(typ))
        print("  K=5 logodds rule types:", dict(typ2))
        results["K5_n_rules"] = int(len(G))
    json.dump(recs, open(f"results/q2_records_K{K}.json", "w"))
json.dump(results, open("results/q2_results.json", "w"), indent=2)
