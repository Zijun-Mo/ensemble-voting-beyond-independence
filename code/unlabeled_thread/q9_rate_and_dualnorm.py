"""Q9 (completion session).
(1) Asymptotic rate of the confidence-set gap  L_pop(maj) - L_conf(maj)  for K=5 at large n
    (n_l = 3200 .. 204800, n_u = 25 n_l), radii as in q7.  Fit log-log slope.
(2) Dual-norm constant C_K(g) = ||lambda*||_1 (marginal multipliers) and ||mu*||_inf (pair multipliers) of the
    sharp-lower-bound LP over E(m,p), for majority and random self-dual LTFs, K = 3..7, many random laws
    (incl. sparse Dirichlet laws).  Conjecture tested:  ||lambda*||_1 <= K  (dictator: 1).
"""
import json
import numpy as np
from scipy.optimize import linprog
from core import *
from q7_finite_sample_and_symK5 import conf_lp, pop_lp  # noqa  (re-uses K=5 globals)

rng = np.random.default_rng(99)
out = {}
K = 5; Zk = patterns(K); FL = flip_index(K); R = reps(K); n = 2 ** K
maj = (Zk.sum(1) >= 3).astype(float)
delta = 0.05
rows = []
for t in range(12):
    alpha, beta, eta = random_joint(K, rng, acc_range=(0.6, 0.85), diff_scale=rng.uniform(0, 3))
    pi = alpha + beta; p = marginals(pi, K)
    if np.any(p <= 0.5):
        continue
    phi = phi_from(alpha, beta, K); m = pair_masses(phi, K)[R]
    Lp, mu, lam = pop_lp(maj, p, pair_masses(phi, K))
    for n_l in (3200, 12800, 51200, 204800):
        n_u = 25 * n_l
        eps_p = np.sqrt(np.log(2 * K / (delta / 2)) / (2 * n_l))
        eps_m = np.sqrt(2 * (len(R) * np.log(2) + np.log(2 / (delta / 2))) / n_u)
        p_hat = rng.binomial(n_l, p) / n_l          # marginal-only labels (independent samples per classifier)
        phi_hat = rng.multinomial(n_u, phi) / n_u
        m_hat = (phi_hat + phi_hat[FL])[R]
        Lc = conf_lp(maj, p_hat, m_hat, eps_p, eps_m)
        rows.append((t, n_l, Lp - Lc, 2 * (eps_p * np.abs(lam).sum() + eps_m * np.abs(mu).max())))
rows = np.array(rows)
ns = np.unique(rows[:, 1])
mg = [rows[rows[:, 1] == v, 2].mean() for v in ns]
slope = np.polyfit(np.log(ns), np.log(mg), 1)[0]
out["rate"] = dict(n_l=ns.tolist(), mean_gap=[float(x) for x in mg], loglog_slope=float(slope),
                   all_valid=bool(np.all(rows[:, 2] >= -1e-9)), all_le_bound=bool(np.all(rows[:, 2] <= rows[:, 3] + 1e-9)))
print(out["rate"], flush=True)

# (2) dual norms
dn = {}
for KK in (3, 4, 5, 6, 7):
    Z = patterns(KK)
    G, _ = enumerate_selfdual_ltfs(KK, n_random=3000) if KK <= 5 else (None, None)
    worst_lam, worst_mu, cnt = 0, 0, 0
    worst_lam_maj = 0
    for t in range(150 if KK <= 5 else 60):
        if t % 2:
            pi = mixture_of_products(KK, rng, diff_scale=rng.uniform(0, 3))
        else:
            pi = rng.dirichlet(np.ones(2 ** KK) * 10 ** rng.uniform(-1.5, 0.3))
            q = rng.uniform(0.55, 0.95, size=KK)
            pi = 0.5 * pi + 0.5 * np.prod(np.where(Z == 1, q, 1 - q), axis=1)
        p = marginals(pi, KK); m = pi + pi[flip_index(KK)]
        A, b = set_even(p, m, KK)
        cands = [(Z.sum(1) > KK / 2).astype(float)] if KK % 2 else []
        if G is not None:
            cands += [G[rng.integers(len(G))].astype(float) for _ in range(5)]
        else:
            for _ in range(5):
                w = rng.exponential(size=KK)
                tab = ltf_table(w, KK)
                if tab is not None:
                    cands.append(tab.astype(float))
        for j, g in enumerate(cands):
            res = linprog(g, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
            y = res.eqlin.marginals
            nl = np.abs(y[2 ** (KK - 1):]).sum(); nm = np.abs(y[:2 ** (KK - 1)]).max()
            worst_lam = max(worst_lam, nl); worst_mu = max(worst_mu, nm); cnt += 1
            if j == 0 and KK % 2:
                worst_lam_maj = max(worst_lam_maj, nl)
    dn[KK] = dict(n_lps=cnt, max_lam_l1=float(worst_lam), max_mu_inf=float(worst_mu), max_lam_l1_majority=float(worst_lam_maj))
    print(KK, dn[KK], flush=True)
out["dual_norms"] = dn
json.dump(out, open("results/q9_results.json", "w"), indent=1)
