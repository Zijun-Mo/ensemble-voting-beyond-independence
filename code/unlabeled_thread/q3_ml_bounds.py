"""Q3: sharp majority-vote error bounds from the first two (weighted) moments of the error count vs
existing PAC-Bayes oracle bounds (first-order 2E[L], tandem 4E[tandem], C-bound = optimised
Chebyshev-Cantelli bound of Wu et al. 2021 Thm 7).

Notation (uniform weights rho_i = 1/K):  W = (1/K) sum_i 1{h_i errs}.
   mu = E_rho[L(h)] = E W,   s = E_{rho^2}[L(h,h')] (diagonal h=h' INCLUDED, as in Masegosa et al.) = E W^2.
   L(MV) = P(W > 1/2) for odd K (no ties).
Sharp bound given (mu, s):  max P(W > 1/2) over laws of W on {0,1/K,...,1} with E W = mu, E W^2 = s.
Continuum version (W in [0,1], ties counted as errors, i.e. P(W >= 1/2)):
   B*(mu,s) = (s - mu^2)/(s - mu + 1/4)   if s <= mu/2        (= C-bound; Cantelli two-point law)
            = 3 mu - 2 s                  if s >= mu/2        (three-point law on {0, 1/2, 1})
Also: information hierarchy on random instances (K=5): C-bound >= sharp(mu,s) >= sharp(p, pairwise) >=
sharp(p, full unlabeled law) >= true error, for uniform and log-odds weights.
"""
import json
import numpy as np
from scipy.optimize import linprog
from core import *

def c_bound(mu, s):
    return np.where(mu < 0.5, (s - mu ** 2) / (s - mu + 0.25), np.nan)

def cc_param(mu, s, nu):
    """Wu et al. 2021 Thm 7: (s - 2 nu mu + nu^2)/(1/2 - nu)^2, nu < 1/2."""
    return (s - 2 * nu * mu + nu ** 2) / (0.5 - nu) ** 2

def sharp_values(vals, err_mask, mu, s):
    """max P(W in err set) over laws on finite support vals with E W = mu, E W^2 = s."""
    A = np.vstack([np.ones_like(vals), vals, vals ** 2])
    res = linprog(-err_mask.astype(float), A_eq=A, b_eq=[1, mu, s], bounds=(0, None), method="highs")
    return -res.fun if res.status == 0 else np.nan

def sharp_K(K, mu, s):
    n = np.arange(K + 1)
    return sharp_values(n / K, n > K / 2, mu, s)

def sharp_cont(mu, s):
    return np.where(s <= mu / 2, (s - mu ** 2) / (s - mu + 0.25), 3 * mu - 2 * s)

def feasible_K(K, mu, s):
    n = np.arange(K + 1)
    A = np.vstack([np.ones(K + 1), n / K, (n / K) ** 2])
    res = linprog(np.zeros(K + 1), A_eq=A, b_eq=[1, mu, s], bounds=(0, None), method="highs")
    return res.status == 0

if __name__ == '__main__':
    out = {}
    # ---------- (a) check continuum closed form against fine-grid LP, and Wu Thm7 min over nu == C-bound
    grid = np.linspace(0, 1, 2001)
    errs = []; errs_cc = []
    rng = np.random.default_rng(3)
    for _ in range(400):
        mu = rng.uniform(0.01, 0.49)
        s = rng.uniform(mu ** 2, mu)
        lp = sharp_values(grid, grid >= 0.5, mu, s)
        errs.append(abs(lp - sharp_cont(mu, s)))
        nu_star = mu - (s - mu ** 2) / (0.5 - mu)          # Wu et al. 2021, optimal offset (can be very negative)
        nus = np.concatenate([np.linspace(-50, 0.4999, 200001), [nu_star]])
        errs_cc.append(abs(np.min(cc_param(mu, s, nus)) - c_bound(mu, s)) + abs(cc_param(mu, s, nu_star) - c_bound(mu, s)))
    out["continuum_closed_form_max_err_vs_gridLP"] = float(np.max(errs))
    out["wu_thm7_min_over_nu_minus_cbound_max_abs"] = float(np.max(errs_cc))
    print(f"continuum closed form vs grid LP: max err {np.max(errs):.2e};  min_nu Thm7 vs C-bound: max |diff| {np.max(errs_cc):.2e}")

    # ---------- (b) grid comparison for K in {3,5,11,21,101, inf}
    def grid_compare(K, n_mu=49, n_r=40):
        rows = []
        for mu in np.linspace(0.01, 0.49, n_mu):
            for frac in np.linspace(0, 1, n_r + 1)[1:-1]:
                s = mu ** 2 + frac * (mu - mu ** 2)
                if K == "inf":
                    sh = float(sharp_cont(mu, s))
                else:
                    if not feasible_K(K, mu, s):
                        continue
                    sh = sharp_K(K, mu, s)
                cb = float(c_bound(mu, s)); fo = 2 * mu; tb = 4 * s
                best_existing = min(cb, fo, tb, 1.0)
                rows.append([mu, s, sh, cb, fo, tb, best_existing])
        return np.array(rows)

    summary = {}
    for K in (3, 5, 11, 21, 101, "inf"):
        R = grid_compare(K)
        sh, cb, fo, tb, be = R[:, 2], R[:, 3], R[:, 4], R[:, 5], R[:, 6]
        sh = np.minimum(sh, 1.0)
        d = dict(n_points=len(R),
                 mean_sharp=float(sh.mean()), mean_cbound=float(np.minimum(cb, 1).mean()), mean_best_existing=float(be.mean()),
                 mean_gap_best_existing_minus_sharp=float((be - sh).mean()), max_gap_best_existing_minus_sharp=float((be - sh).max()),
                 mean_ratio_sharp_over_best_existing=float(np.mean(sh / be)), min_ratio_sharp_over_best_existing=float(np.min(sh / be)),
                 frac_best_existing_is_sharp=float(np.mean(be - sh < 1e-6)),
                 mean_gap_cbound_minus_sharp=float((np.minimum(cb, 1) - sh).mean()),
                 frac_cbound_is_sharp=float(np.mean(np.minimum(cb, 1) - sh < 1e-6)),
                 frac_region_s_le_mu_over_2=float(np.mean(R[:, 1] <= R[:, 0] / 2)))
        summary[str(K)] = d
        np.save(f"results/q3_grid_K{K}.npy", R)
        print(f"K={K}: n={len(R)} mean sharp {d['mean_sharp']:.4f} | mean best existing {d['mean_best_existing']:.4f} | "
              f"mean gap {d['mean_gap_best_existing_minus_sharp']:.4f} (max {d['max_gap_best_existing_minus_sharp']:.4f}) | "
              f"mean ratio sharp/best {d['mean_ratio_sharp_over_best_existing']:.3f} (min {d['min_ratio_sharp_over_best_existing']:.3f}) | "
              f"best-existing sharp at {d['frac_best_existing_is_sharp']:.3f}; C-bound sharp at {d['frac_cbound_is_sharp']:.3f}")
    out["grid"] = summary

    # A few illustrative points (mu=0.2)
    ill = []
    for K in (3, 5, 11, 21, "inf"):
        for s in (0.05, 0.08, 0.10, 0.14):
            mu = 0.2
            sh = float(sharp_cont(mu, s)) if K == "inf" else (sharp_K(K, mu, s) if feasible_K(K, mu, s) else np.nan)
            ill.append(dict(K=str(K), mu=mu, s=s, sharp=sh, cbound=float(c_bound(mu, s)), first_order=2 * mu, tandem=4 * s))
    out["illustrative_mu0.2"] = ill
    for r in ill:
        print(f"  K={r['K']:>4} mu=0.2 s={r['s']:.2f}: sharp {r['sharp']:.4f}  C {r['cbound']:.4f}  2mu {r['first_order']:.3f}  4s {r['tandem']:.3f}")

    # ---------- (c) information hierarchy on random instances (K=5)
    K = 5
    Z = patterns(K)
    E = 1 - Z  # error indicators
    hier = {"uniform": [], "logodds": []}
    for it in range(400):
        pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3), acc_range=(0.55, 0.9))
        p = marginals(pi, K)
        if np.any(p <= 0.5):
            continue
        r = pairwise_joint(pi, K)
        m = pi + pi[flip_index(K)]
        A_p, b_p = set_pairwise(p, r, K)
        A_e, b_e = set_even(p, m, K)
        for wname, w in (("uniform", np.ones(K)), ("logodds", np.log(p / (1 - p)))):
            rho = w / w.sum()
            Wz = E @ rho
            err = (Wz > 0.5).astype(float)
            mu = float(Wz @ pi); s = float((Wz ** 2) @ pi)
            true = float(err @ pi)
            vals, inv = np.unique(np.round(Wz, 12), return_inverse=True)
            sh_ms = sharp_values(vals, vals > 0.5, mu, s)
            sh_pr = -lp_min(-err, A_p, b_p)
            sh_full = -lp_min(-err, A_e, b_e)
            hier[wname].append([true, sh_full, sh_pr, sh_ms, float(c_bound(mu, s)), 2 * mu, 4 * s])
    hsum = {}
    for wname, rows in hier.items():
        R = np.array(rows)
        names = ["true", "sharp_p_fullunlabeled", "sharp_p_pairwise", "sharp_mu_s", "C_bound", "first_order_2mu", "tandem_4s"]
        hsum[wname] = {nm: float(np.mean(np.minimum(R[:, i], 1))) for i, nm in enumerate(names)}
        hsum[wname]["n"] = len(R)
        hsum[wname]["monotone_chain_holds"] = bool(np.all(np.diff(R[:, :5], axis=1) >= -1e-8))
        hsum[wname]["frac_sharp_mu_s_eq_C"] = float(np.mean(np.abs(R[:, 3] - R[:, 4]) < 1e-6))
        hsum[wname]["median_ratio_C_over_sharp_pairwise"] = float(np.median(R[:, 4] / R[:, 2]))
        np.save(f"results/q3_hierarchy_{wname}.npy", R)
        print(f"hierarchy K=5 [{wname}] n={len(R)}: " + ", ".join(f"{nm} {hsum[wname][nm]:.4f}" for nm in names)
              + f"; chain holds: {hsum[wname]['monotone_chain_holds']}; sharp(mu,s)==C in {hsum[wname]['frac_sharp_mu_s_eq_C']:.3f}")
    out["hierarchy_K5"] = hsum

    # ---------- (d) K=3 closed form of sharp error bound given individual errors e_i and pairwise tandem/disagreement
    K = 3
    Z = patterns(K)
    maxerr = 0
    for _ in range(2000):
        pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3), acc_range=(0.5, 0.95))
        p = marginals(pi, K); e = 1 - p
        r = pairwise_joint(pi, K)
        tand = {(i, j): 1 - p[i] - p[j] + r[(i, j)] for (i, j) in r}
        dis = {(i, j): e[i] + e[j] - 2 * tand[(i, j)] for (i, j) in r}
        T = sum(tand.values())
        d = np.array([(dis[tuple(sorted((i, j)))] + dis[tuple(sorted((i, k)))] - dis[tuple(sorted((j, k)))]) / 2
                      for i, j, k in ((0, 1, 2), (1, 0, 2), (2, 0, 1))])
        cf = min(T, np.min(e + d))
        A_p, b_p = set_pairwise(p, r, K)
        err = (Z.sum(1) <= 1).astype(float)
        lp = -lp_min(-err, A_p, b_p)
        maxerr = max(maxerr, abs(cf - lp))
    out["K3_closed_form_sharp_error_max_err"] = float(maxerr)
    print(f"K=3 sharp MV error = min(sum_i<j tandem_ij, min_i(e_i + d_i)) vs LP: max err {maxerr:.2e}")
    json.dump(out, open("results/q3_results.json", "w"), indent=2)
