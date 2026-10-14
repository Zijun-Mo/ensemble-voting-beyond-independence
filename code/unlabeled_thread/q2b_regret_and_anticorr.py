"""Q2b: (i) negatively-dependent ('disjoint error regions') family where max-min can be MORE aggressive
than log-odds; (ii) randomized minimax-REGRET mixtures of self-dual LTFs as a less conservative robust rule;
(iii) which rule types max-min picks when the FULL unlabeled law (E(m,p)) is used for K=5.
"""
import json, time
from collections import Counter
import numpy as np
from scipy.optimize import linprog
from core import *


rng = np.random.default_rng(77)


def anti_instance(K, rng):
    """U ~ Unif[0,1); classifier i errs iff U in an arc of length e_i starting at a random offset.
    Offsets spread out -> errors (nearly) disjoint; clustered -> overlapping. Exact pi by discretisation."""
    e = rng.uniform(0.1, 0.4, size=K)
    spread = rng.uniform(0, 1)
    starts = (np.arange(K) / K * spread + rng.uniform(0, 1)) % 1.0
    N = 4000
    u = (np.arange(N) + 0.5) / N
    err = np.zeros((N, K), bool)
    for i in range(K):
        err[:, i] = ((u - starts[i]) % 1.0) < e[i]
    z = (~err).astype(int)
    idx = (z * (2 ** np.arange(K))).sum(1)
    pi = np.bincount(idx, minlength=2 ** K) / N
    return pi


def minimax_regret_mixture(G, A, b):
    nr, n = G.shape
    nd = A.shape[0]
    # vars: t, omega (nr), y (nr*nd)
    nv = 1 + nr + nr * nd
    c = np.zeros(nv); c[0] = 1
    A_ub, b_ub = [], []
    for h in range(nr):
        row = np.zeros(nv); row[0] = -1; row[1 + nr + h * nd: 1 + nr + (h + 1) * nd] = b
        A_ub.append(row); b_ub.append(0)
    blocks = []
    for h in range(nr):
        # -(A^T y_h) - G^T omega <= -G_h
        M = np.zeros((n, nv))
        M[:, 1:1 + nr] = -G.T
        M[:, 1 + nr + h * nd: 1 + nr + (h + 1) * nd] = -A.T
        blocks.append(M)
        b_ub.extend(list(-G[h]))
    A_ub = np.vstack([np.array(A_ub)] + blocks)
    A_eq = np.zeros((1, nv)); A_eq[0, 1:1 + nr] = 1
    bounds = [(None, None)] + [(0, None)] * nr + [(None, None)] * (nr * nd)
    res = linprog(c, A_ub=A_ub, b_ub=np.array(b_ub), A_eq=A_eq, b_eq=[1.0], bounds=bounds, method="highs")
    if res.status != 0:
        return None, None
    return res.fun, res.x[1:1 + nr]


def det_regret(G, A, b, cand):
    """max regret of each candidate deterministic rule k vs all rules h (LPs)."""
    out = {}
    for k in cand:
        worst = 0.0
        for h in range(len(G)):
            if h == k:
                continue
            res = linprog(-(G[h] - G[k]).astype(float), A_eq=A, b_eq=b, bounds=(0, None), method="highs")
            worst = max(worst, -res.fun)
        out[k] = worst
    return out


if __name__ == '__main__':
    results = {}
    for K in (3, 5):
        G, W = enumerate_selfdual_ltfs(K, n_random=50000 if K == 5 else 2000)
        Gf = G.astype(float)
        Z = patterns(K)
        nr = len(G)
        table_index = {tuple(G[k]): k for k in range(nr)}
        maj = table_index[tuple((Z.sum(1) > K / 2).astype(int))]
        dicts = [table_index[tuple(Z[:, i].astype(int))] for i in range(K)]
        types = {k: tuple(sorted(canonical_weights(G[k], K, maxw=5), reverse=True)) for k in range(nr)}
        for fam in ("independent", "mixture", "bloc", "anti"):
            t0 = time.time()
            recs = []
            n_inst = 150 if K == 5 else 300
            while len(recs) < n_inst:
                pi = anti_instance(K, rng) if fam == "anti" else family_instance(K, fam, rng)
                p = marginals(pi, K)
                if np.any(p <= 0.5):
                    continue
                r = pairwise_joint(pi, K)
                m = pi + pi[flip_index(K)]
                A_p, b_p = set_pairwise(p, r, K)
                A_e, b_e = set_even(p, m, K)
                true = Gf @ pi
                lo_pw = np.array([lp_min(Gf[k], A_p, b_p) for k in range(nr)])
                lo_full = np.array([lp_min(Gf[k], A_e, b_e) for k in range(nr)])
                k_lo = table_index[tuple(ltf_table(np.log(p / (1 - p)), K))]
                k_rpw, k_rfull = int(np.argmax(lo_pw)), int(np.argmax(lo_full))
                reg_val, omega = minimax_regret_mixture(Gf, A_p, b_p)
                reg_true = float(omega @ true)
                k_reg_det = int(np.argmax(omega))
                rec = dict(true_lo=true[k_lo], true_rpw=true[k_rpw], true_rfull=true[k_rfull], true_maj=true[maj],
                           true_best=max(true[d] for d in dicts), true_reg=reg_true, true_reg_det=true[k_reg_det],
                           oracle=true.max(), lo_pw_lo=lo_pw[k_lo], lo_pw_rpw=lo_pw[k_rpw], lo_pw_maj=lo_pw[maj],
                           lo_full_lo=lo_full[k_lo], lo_full_rfull=lo_full[k_rfull], minimax_regret=reg_val,
                           type_rpw=types[k_rpw], type_rfull=types[k_rfull], type_lo=types[k_lo],
                           type_regdet=types[k_reg_det], differ=k_rpw != k_lo,
                           rob_maj_lo_dict=(k_rpw == maj and k_lo in dicts))
                recs.append(rec)
            R = recs
            s = {key: float(np.mean([r[key] for r in R])) for key in
                 ("true_lo", "true_rpw", "true_rfull", "true_maj", "true_best", "true_reg", "true_reg_det", "oracle",
                  "lo_pw_lo", "lo_pw_rpw", "lo_pw_maj", "lo_full_lo", "lo_full_rfull", "minimax_regret")}
            s["n"] = len(R)
            s["frac_differ_rpw_lo"] = float(np.mean([r["differ"] for r in R]))
            s["frac_robust_majority_while_logodds_dictator"] = float(np.mean([r["rob_maj_lo_dict"] for r in R]))
            s["types_rpw"] = {str(k): v for k, v in Counter(r["type_rpw"] for r in R).items()}
            s["types_rfull"] = {str(k): v for k, v in Counter(r["type_rfull"] for r in R).items()}
            s["types_regret_det"] = {str(k): v for k, v in Counter(r["type_regdet"] for r in R).items()}
            s["true_reg_minus_lo"] = float(np.mean([r["true_reg"] - r["true_lo"] for r in R]))
            s["true_rpw_minus_lo"] = float(np.mean([r["true_rpw"] - r["true_lo"] for r in R]))
            s["frac_true_rpw_gt_lo"] = float(np.mean([r["true_rpw"] > r["true_lo"] + 1e-12 for r in R]))
            s["frac_true_rpw_lt_lo"] = float(np.mean([r["true_rpw"] < r["true_lo"] - 1e-12 for r in R]))
            results[f"K{K}_{fam}"] = s
            print(f"K={K} {fam:12s} n={len(R)} ({time.time()-t0:.0f}s): true acc  logodds {s['true_lo']:.4f} | maxmin_pw {s['true_rpw']:.4f} "
                  f"| maxmin_full {s['true_rfull']:.4f} | regret-mix {s['true_reg']:.4f} | regret-det {s['true_reg_det']:.4f} | maj {s['true_maj']:.4f} "
                  f"| best {s['true_best']:.4f} | oracle {s['oracle']:.4f}")
            print(f"     worst(pw): logodds {s['lo_pw_lo']:.4f} maxmin {s['lo_pw_rpw']:.4f} maj {s['lo_pw_maj']:.4f}; worst(full): logodds {s['lo_full_lo']:.4f} maxmin {s['lo_full_rfull']:.4f}; "
                  f"minimax regret {s['minimax_regret']:.4f}; differ {s['frac_differ_rpw_lo']:.3f}; rob=maj & lo=dict {s['frac_robust_majority_while_logodds_dictator']:.3f}; "
                  f"P(rpw>lo) {s['frac_true_rpw_gt_lo']:.3f} P(rpw<lo) {s['frac_true_rpw_lt_lo']:.3f}")
            print(f"     types maxmin_pw {s['types_rpw']}\n     types maxmin_full {s['types_rfull']}\n     types regret_det {s['types_regret_det']}")
    json.dump(results, open("results/q2b_results.json", "w"), indent=2)
