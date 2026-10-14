"""Q1: what is identified about pi = law(Z) from (law of f, marginals p)?

Checks (all on random joint laws of (Y,Z) with class-dependent error patterns):
 A. Prior unknown: identified set of pi == E(m,p) = {pi>=0: pi(z)+pi(1-z)=m(z), marginals p}.
    Test: range of random linear functionals c.pi over direct (alpha,beta) LP == over E(m,p).
 B. Prior known (eta): identified set == E(m,p) cap {TV(pi,phi)<=1-eta} cap {TV(pi,phi~)<=eta}.
 C. K=3: E(m,p) == P(p,r) (pairwise set) and the closed form of the majority interval in observables
    d_i = P(f_i is the lone dissenter): [max{p_i-d_i, sum p-2+D}, min{p_i+d_i, sum p - D}].
 D. Third-order t = P(Z1=Z2=Z3=1): interval width with prior unknown vs known.
 E. Z independent of Y, eta != 1/2: pi recovered exactly as (phi - phi~)/(2 eta - 1) + even part.
 F. Strong duality: sharp lower bound == max_lambda 1/2[1 + lam.(2p-1) - E_phi|1 - G lam.s|]
    and the one-dimensional 'margin certificate' along lam = c*w for weighted votes.
 G. K=5: width of the majority / weighted-vote intervals under pairwise-only vs full-unlabeled vs +prior.
"""
import json, sys, time
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from core import *

rng = np.random.default_rng(12345)
out = {}

# ---------------- A & B ----------------
def check_AB(K, n_inst=150, n_obj=4):
    maxdiff_A, maxdiff_B = 0.0, 0.0
    shrink = []
    for _ in range(n_inst):
        alpha, beta, eta = random_joint(K, rng, diff_scale=rng.uniform(0, 3))
        pi = alpha + beta
        phi = phi_from(alpha, beta, K)
        p = marginals(pi, K)
        m = pair_masses(phi, K)
        A_e, b_e = set_even(p, m, K)
        for _ in range(n_obj):
            c = rng.normal(size=2 ** K)
            r1 = direct_ab_range(c, phi, p, K)
            r2 = lp_range(c, A_e, b_e)
            maxdiff_A = max(maxdiff_A, abs(r1[0] - r2[0]), abs(r1[1] - r2[1]))
            r3 = direct_ab_range(c, phi, p, K, eta=eta)
            r4 = tv_range(c, phi, p, K, eta)
            maxdiff_B = max(maxdiff_B, abs(r3[0] - r4[0]), abs(r3[1] - r4[1]))
            shrink.append((r3[1] - r3[0]) / max(r1[1] - r1[0], 1e-12))
    return maxdiff_A, maxdiff_B, float(np.mean(shrink)), float(np.median(shrink))

for K in (3, 4, 5):
    t0 = time.time()
    a, b, ms, md = check_AB(K, n_inst=120 if K < 5 else 60)
    out[f"AB_K{K}"] = dict(maxdiff_prior_unknown=a, maxdiff_prior_known=b,
                          mean_width_ratio_known_over_unknown=ms, median_width_ratio=md)
    print(f"[A/B] K={K}: max|direct - E(m,p)| = {a:.2e}; max|direct(eta) - TV-set| = {b:.2e}; "
          f"width ratio (eta known / unknown): mean {ms:.3f}, median {md:.3f}  ({time.time()-t0:.1f}s)")

# ---------------- C: K=3 closed form ----------------
K = 3
Z = patterns(K)
maj = (Z.sum(1) >= 2).astype(float)
maxerr_lo = maxerr_hi = maxerr_set = 0.0
n_c = 3000
for _ in range(n_c):
    pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3), acc_range=(0.5, 0.95))
    p = marginals(pi, K)
    r = pairwise_joint(pi, K)
    m = pi + pi[flip_index(K)]
    # lone dissent probabilities: pattern pairs where i differs from the other two
    d = np.array([sum(pi[v] for v in range(8) if (Z[v, i] != Z[v, (i + 1) % 3]) and Z[v, (i + 1) % 3] == Z[v, (i + 2) % 3]) for i in range(3)])
    D = d.sum()
    lo_cf = max(np.max(p - d), p.sum() - 2 + D)
    hi_cf = min(np.min(p + d), p.sum() - D)
    A_e, b_e = set_even(p, m, K)
    lo1, hi1 = lp_range(maj, A_e, b_e)
    A_p, b_p = set_pairwise(p, r, K)
    lo2, hi2 = lp_range(maj, A_p, b_p)
    maxerr_lo = max(maxerr_lo, abs(lo_cf - lo1)); maxerr_hi = max(maxerr_hi, abs(hi_cf - hi1))
    maxerr_set = max(maxerr_set, abs(lo1 - lo2), abs(hi1 - hi2))
out["C_K3_closed_form"] = dict(n=n_c, max_abs_err_lower=maxerr_lo, max_abs_err_upper=maxerr_hi,
                               max_abs_diff_even_vs_pairwise=maxerr_set)
print(f"[C] K=3 closed form vs LP over {n_c} instances: lower err {maxerr_lo:.2e}, upper err {maxerr_hi:.2e}; "
      f"E(m,p) vs P(p,r): {maxerr_set:.2e}")

# ---------------- D: third-order t ----------------
K = 3
t_obj = np.prod(patterns(K), axis=1)
widths_unknown, widths_known, n_strict = [], [], 0
for _ in range(500):
    alpha, beta, eta = random_joint(K, rng, diff_scale=rng.uniform(0, 3))
    pi = alpha + beta; phi = phi_from(alpha, beta, K); p = marginals(pi, K)
    w1 = direct_ab_range(t_obj, phi, p, K)
    w2 = direct_ab_range(t_obj, phi, p, K, eta=eta)
    widths_unknown.append(w1[1] - w1[0]); widths_known.append(w2[1] - w2[0])
    n_strict += (w2[1] - w2[0]) < (w1[1] - w1[0]) - 1e-9
widths_unknown = np.array(widths_unknown); widths_known = np.array(widths_known)
out["D_t_identification"] = dict(n=500, mean_width_prior_unknown=float(widths_unknown.mean()),
                                 mean_width_prior_known=float(widths_known.mean()),
                                 frac_positive_width_prior_unknown=float(np.mean(widths_unknown > 1e-9)),
                                 frac_positive_width_prior_known=float(np.mean(widths_known > 1e-9)),
                                 frac_strictly_shrunk_by_prior=n_strict / 500)
print(f"[D] t=P(Z1=Z2=Z3=1): mean width prior unknown {widths_unknown.mean():.4f}, prior known {widths_known.mean():.4f}; "
      f"positive width in {np.mean(widths_known > 1e-9):.3f} of instances even with prior known; "
      f"prior shrinks interval in {n_strict/500:.3f}")

# ---------------- E: Z independent of Y ----------------
errs = []
for K in (3, 5):
    for _ in range(200):
        pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3))
        eta = rng.uniform(0.2, 0.8)
        alpha, beta = eta * pi, (1 - eta) * pi
        phi = phi_from(alpha, beta, K)
        p = marginals(pi, K)
        fl = flip_index(K)
        # eta from one classifier: 2P(f_i=1)-1 = (2eta-1)(2p_i-1)
        i = int(np.argmax(np.abs(2 * p - 1)))
        Pf1 = patterns(K)[:, i] @ phi
        eta_hat = 0.5 * (1 + (2 * Pf1 - 1) / (2 * p[i] - 1))
        odd = (phi - phi[fl]) / (2 * eta_hat - 1)
        pi_hat = 0.5 * ((phi + phi[fl]) + odd)
        errs.append(np.max(np.abs(pi_hat - pi)))
out["E_independent_case_recovery_maxerr"] = float(np.max(errs))
print(f"[E] Z indep of Y, eta != 1/2: max |pi_hat - pi| = {np.max(errs):.2e}")

# ---------------- F: duality / margin certificate ----------------
def dual_lp(g, phi, p, K):
    """max_lam 1/2[1 + lam.(2p-1) - sum_v phi_v |1 - G_v lam.s_v|] as LP."""
    n = 2 ** K
    S = 2 * patterns(K) - 1
    Gs = 2 * g - 1
    # vars lam (K free), u (n >= |1 - G lam.s|)
    c = np.concatenate([-(2 * p - 1), phi])
    M = Gs[:, None] * S
    A_ub = np.vstack([np.hstack([-M, -np.eye(n)]), np.hstack([M, -np.eye(n)])])
    b_ub = np.concatenate([-np.ones(n), np.ones(n)])
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(None, None)] * K + [(0, None)] * n, method="highs")
    return 0.5 * (1 - res.fun), res.x[:K]

def ray_certificate(w, phi, p, K):
    S = 2 * patterns(K) - 1
    marg = S @ w
    g = (marg > 0).astype(float)
    f = lambda c: -dual_value(c * w, g, phi, p, K)
    res = minimize_scalar(f, bounds=(0, 50 / np.abs(w).sum()), method="bounded", options=dict(xatol=1e-10))
    return -res.fun

dual_err = 0.0
ray_gap = []
for K in (3, 5):
    G5, W5 = enumerate_selfdual_ltfs(K, n_random=20000)
    for _ in range(150):
        alpha, beta, eta = random_joint(K, rng, diff_scale=rng.uniform(0, 3))
        pi = alpha + beta; phi = phi_from(alpha, beta, K); p = marginals(pi, K)
        A_e, b_e = set_even(p, pair_masses(phi, K), K)
        k = rng.integers(len(G5))
        g = G5[k].astype(float)
        lo, hi = lp_range(g, A_e, b_e)
        dv, lam = dual_lp(g, phi, p, K)
        dual_err = max(dual_err, abs(dv - lo))
        # ray certificate along the rule's own weights (log-odds weights here)
        w = np.log(np.clip(p, 0.501, 0.999) / (1 - np.clip(p, 0.501, 0.999)))
        gw = ltf_table(w, K)
        lo_w, _ = lp_range(gw.astype(float), A_e, b_e)
        ray_gap.append(lo_w - ray_certificate(w, phi, p, K))
ray_gap = np.array(ray_gap)
out["F_duality"] = dict(max_abs_err_dual_vs_primal=dual_err,
                        ray_certificate_gap_mean=float(ray_gap.mean()), ray_certificate_gap_median=float(np.median(ray_gap)),
                        ray_certificate_gap_max=float(ray_gap.max()), frac_ray_exact=float(np.mean(ray_gap < 1e-6)))
print(f"[F] dual formula vs LP: max err {dual_err:.2e}; ray certificate (lam = c*w) gap to sharp: "
      f"mean {ray_gap.mean():.4f}, median {np.median(ray_gap):.4f}, max {ray_gap.max():.4f}, exact in {np.mean(ray_gap<1e-6):.3f}")

# ---------------- G: K=5 widths ----------------
K = 5
Z = patterns(K)
maj = (Z.sum(1) >= 3).astype(float)
rows = []
for _ in range(300):
    alpha, beta, eta = random_joint(K, rng, diff_scale=rng.uniform(0, 3))
    pi = alpha + beta; phi = phi_from(alpha, beta, K); p = marginals(pi, K)
    r = pairwise_joint(pi, K)
    A_p, b_p = set_pairwise(p, r, K)
    A_e, b_e = set_even(p, pair_masses(phi, K), K)
    a1 = lp_range(maj, A_p, b_p); a2 = lp_range(maj, A_e, b_e); a3 = tv_range(maj, phi, p, K, eta)
    truth = maj @ pi
    rows.append([a1[0], a1[1], a2[0], a2[1], a3[0], a3[1], truth])
rows = np.array(rows)
wp, we, wt = rows[:, 1] - rows[:, 0], rows[:, 3] - rows[:, 2], rows[:, 5] - rows[:, 4]
out["G_K5_majority_interval_width"] = dict(n=len(rows), pairwise_mean=float(wp.mean()), full_unlabeled_mean=float(we.mean()),
                                          full_plus_prior_mean=float(wt.mean()),
                                          ratio_full_over_pairwise_mean=float(np.mean(we / wp)),
                                          lower_gain_full_minus_pairwise_mean=float(np.mean(rows[:, 2] - rows[:, 0])),
                                          truth_inside_all=bool(np.all((rows[:, 4] <= rows[:, 6] + 1e-9) & (rows[:, 6] <= rows[:, 5] + 1e-9))))
print(f"[G] K=5 majority interval widths: pairwise {wp.mean():.4f}, full unlabeled {we.mean():.4f}, +prior {wt.mean():.4f}; "
      f"ratio full/pairwise mean {np.mean(we/wp):.3f}; lower-bound gain {np.mean(rows[:,2]-rows[:,0]):.4f}")
np.save("results/q1_K5_intervals.npy", rows)
json.dump(out, open("results/q1_results.json", "w"), indent=2)
