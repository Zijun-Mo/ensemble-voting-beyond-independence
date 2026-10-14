"""Q7 (completion session).
(A) Finite-sample theorem check.  Confidence set
      C = { pi in simplex : |Z^T pi - p_hat|_inf <= eps_p ,  sum_{v in reps} |pi(v)+pi(vbar) - m_hat(v)| <= eps_m }
    with eps_p = sqrt(log(2K/d1)/(2 n_l))  (Hoeffding + union),  eps_m = Weissman L1 radius
    sqrt(2 (ln(2^{2^{K-1}} - 2) + ln(1/d2)) / n_u)  for the 2^{K-1}-cell law m.
    Claim (proved in notes):  on the event {(p,m) in radii},
        L(true) - 2 (eps_p ||lam*||_1 + eps_m ||mu*||_inf)  <=  L_conf  <=  L(true)  <=  acc(g),
    where (mu*, lam*) is any optimal dual of the population LP  min{g.pi : pi in E(m,p)}.
    Also: max-min rule chosen on C has true worst case >= max_g L(g) - max_g gap_g.
    Check coverage, the gap inequality, and 1/sqrt(n) scaling.
    Also report observed dual norms ||lam*||_1, ||mu*||_inf per K=5 rule type (constant C_K(g)).
(B) Symmetric K=5 closed forms: for p_i=p, r_ij=r, the sharp majority lower bound over P(p,r) equals the
    exchangeable (number-of-errors) two-moment closed form (Thm 6 of the earlier notes) -> check vs LP;
    max-min value = max{p, A3, A5} with A3 = K=3 closed form, A5 = Thm-6 form.
"""
import json
import numpy as np
from math import comb
from scipy.optimize import linprog
from core import *

rng = np.random.default_rng(4242)
out = {}
K = 5
Zk = patterns(K); FL = flip_index(K); R = reps(K); n = 2 ** K
G, W = enumerate_selfdual_ltfs(K, n_random=20000); Gf = G.astype(float)
types = [canonical_weights(G[k], K) for k in range(len(G))]


def pop_lp(g, p, m):
    A, b = set_even(p, m, K)
    res = linprog(g, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    y = res.eqlin.marginals
    return res.fun, y[:len(R)], y[len(R):]   # value, mu (pair rows), lam (marginals)


def conf_lp(g, p_hat, m_hat, eps_p, eps_m):
    # vars: pi (n), u (len R) ;  |pair(pi) - m_hat| <= u, sum u <= eps_m
    nr = len(R)
    P = np.zeros((nr, n))
    for k, v in enumerate(R):
        P[k, v] += 1; P[k, FL[v]] += 1
    A_ub = np.vstack([np.hstack([Zk.T, np.zeros((K, nr))]), np.hstack([-Zk.T, np.zeros((K, nr))]),
                      np.hstack([P, -np.eye(nr)]), np.hstack([-P, -np.eye(nr)]),
                      np.concatenate([np.zeros(n), np.ones(nr)])[None, :]])
    b_ub = np.concatenate([p_hat + eps_p, -(p_hat - eps_p), m_hat, -m_hat, [eps_m]])
    A_eq = np.concatenate([np.ones(n), np.zeros(nr)])[None, :]
    res = linprog(np.concatenate([g, np.zeros(nr)]), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0],
                  bounds=(0, None), method="highs")
    return res.fun if res.status == 0 else np.nan


if __name__ == "__main__":
    delta = 0.05
    recs = []
    for t in range(40):
        alpha, beta, eta = random_joint(K, rng, acc_range=(0.6, 0.85), diff_scale=rng.uniform(0, 3))
        pi = alpha + beta; p = marginals(pi, K)
        if np.any(p <= 0.5):
            continue
        phi = phi_from(alpha, beta, K); m = pair_masses(phi, K)[R]
        pops = [pop_lp(g, p, pair_masses(phi, K)) for g in Gf]
        Lpop = np.array([v[0] for v in pops])
        for n_l, n_u in ((200, 5000), (800, 20000), (3200, 80000)):
            eps_p = np.sqrt(np.log(2 * K / (delta / 2)) / (2 * n_l))
            eps_m = np.sqrt(2 * (len(R) * np.log(2) + np.log(2 / (delta / 2))) / n_u)
            # sample
            idx_l = rng.choice(n, size=n_l, p=pi)          # correctness patterns of labeled points
            p_hat = Zk[idx_l].mean(0)
            sgn = rng.random(n_u) < eta
            a_n = alpha / alpha.sum(); b_n = beta / beta.sum()
            zz = np.where(sgn, rng.choice(n, size=n_u, p=a_n), FL[rng.choice(n, size=n_u, p=b_n)])  # f pattern
            phi_hat = np.bincount(zz, minlength=n) / n_u
            m_hat = (phi_hat + phi_hat[FL])[R]
            event = (np.abs(p_hat - p).max() <= eps_p) and (np.abs(m_hat - m).sum() <= eps_m)
            Lc = np.array([conf_lp(g, p_hat, m_hat, eps_p, eps_m) for g in Gf])
            gaps = Lpop - Lc
            bnd = np.array([2 * (eps_p * np.abs(v[2]).sum() + eps_m * np.abs(v[1]).max()) for v in pops])
            k_hat = int(np.nanargmax(Lc))
            recs.append(dict(n_l=n_l, event=bool(event), valid=bool(np.all(Lc <= Lpop + 1e-9)),
                             gap_le_bound=bool(np.all(gaps <= bnd + 1e-9)), mean_gap=float(np.mean(gaps)),
                             maj_gap=float(gaps[[i for i, w in enumerate(types) if w == (1, 1, 1, 1, 1)][0]]),
                             mean_bound=float(np.mean(bnd)), eps_p=eps_p, eps_m=eps_m,
                             maxmin_regret=float(Lpop.max() - Lpop[k_hat]),
                             lam_norm=[float(np.abs(v[2]).sum()) for v in pops], mu_norm=[float(np.abs(v[1]).max()) for v in pops]))
    summ = {}
    for n_l in (200, 800, 3200):
        rr = [r for r in recs if r["n_l"] == n_l]
        summ[n_l] = dict(n=len(rr), frac_event=np.mean([r["event"] for r in rr]), frac_valid=np.mean([r["valid"] for r in rr]),
                         frac_gap_le_bound=np.mean([r["gap_le_bound"] for r in rr]), mean_gap=np.mean([r["mean_gap"] for r in rr]),
                         mean_maj_gap=np.mean([r["maj_gap"] for r in rr]), mean_bound=np.mean([r["mean_bound"] for r in rr]),
                         eps_p=rr[0]["eps_p"], eps_m=rr[0]["eps_m"], mean_maxmin_regret=np.mean([r["maxmin_regret"] for r in rr]),
                         max_maxmin_regret=np.max([r["maxmin_regret"] for r in rr]))
        print("A", n_l, {k: round(float(v), 4) for k, v in summ[n_l].items()}, flush=True)
    out["A_finite_sample"] = {str(k): {kk: float(vv) for kk, vv in v.items()} for k, v in summ.items()}
    lam = np.array([r["lam_norm"] for r in recs]); mu = np.array([r["mu_norm"] for r in recs])
    bytype = {}
    for k, w in enumerate(types):
        key = str(tuple(sorted(w, reverse=True)))
        bytype.setdefault(key, []).append(k)
    out["A_dual_norms_by_type"] = {key: dict(lam_l1_max=float(lam[:, ks].max()), lam_l1_mean=float(lam[:, ks].mean()),
                                             mu_inf_max=float(mu[:, ks].max())) for key, ks in bytype.items()}
    print(out["A_dual_norms_by_type"])

    # ------------------------------------------------------------------ (B) symmetric K=5 closed forms
    def thm6(Kk, mu_, s_):
        r = (Kk + 1) // 2; m1, m2 = Kk * mu_, Kk * Kk * s_
        c = [1.0, ((r + Kk) * m1 - m2) / (r * Kk)]
        for j in range(0, r - 1):
            c.append((m2 - (2 * j + 1) * m1 + j * (j + 1)) / ((r - j) * (r - j - 1)))
        return min(c)


    def A3(p, rr):
        e = 1 - p; T = 1 - 2 * p + rr; D = 2 * (e - T)     # D = P(disagree) for a pair
        return 1 - min(3 * T, e + D / 2)                      # K=3 closed form, symmetric (d_i = D/2)


    maj5 = (Zk.sum(1) >= 3).astype(float); maj3 = (Zk[:, :3].sum(1) >= 2).astype(float)
    errs5, errs3, errmm, nn = [], [], [], 0
    for p in np.linspace(0.52, 0.98, 24):
        for rho in np.linspace(-0.3, 1.0, 40):
            rr = p * p + rho * p * (1 - p)
            A, b = set_pairwise(np.full(K, p), {(i, j): rr for i in range(K) for j in range(i + 1, K)}, K)
            v5 = lp_min(maj5, A, b)
            if v5 is None:
                continue
            v3 = lp_min(maj3, A, b)
            e = 1 - p; T = 1 - 2 * p + rr
            mu_, s_ = e, (K * e + K * (K - 1) * T) / K ** 2
            c5 = 1 - thm6(K, mu_, s_)
            errs5.append(abs(v5 - c5)); errs3.append(abs(v3 - A3(p, rr)))
            lows = [lp_min(g, A, b) for g in Gf]
            errmm.append(abs(max(lows) - max(p, A3(p, rr), c5))); nn += 1
    out["B_symmetric_K5"] = dict(n=nn, maj5_closed_vs_lp_maxerr=float(max(errs5)), maj3_closed_vs_lp_maxerr=float(max(errs3)),
                                 maxmin_value_closed_vs_lp_maxerr=float(max(errmm)))
    print("B", out["B_symmetric_K5"])
    json.dump(out, open("results/q7_results.json", "w"), indent=1, default=float)
