"""Q3b / Q2c: closed forms and structural checks.

1. Finite-K closed form of the sharp MV-error bound given (mu, s) [uniform weights, odd K, r=(K+1)/2, nu=K W]:
     U = min{1, E[nu (r+K-nu)]/(rK), min_{j=0..r-2} E[(nu-j)(nu-j-1)]/((r-j)(r-j-1))}
   (dual quadratic polynomials touching {0,r,K} or {j,j+1,r}); check vs LP.
   Continuum: U = min{C-bound, 3mu-2s} -> check  'L(MV) <= 3 mu - 2 s = E[L] + E[disagreement]' is always valid.
2. Weighted vote: sharp bound given weighted (mu_rho, s_rho) over {0,1}^K equals the LP over the value set
   V_rho = {sum_{i in S} rho_i}; generalised closed form with V_rho; check vs 2^K LP for K=3,5.
3. K=3 identities: D = P(not unanimous) = sum e_i - sum_{i<j} T_ij (T = tandem loss); max-min picks majority iff
   sum_{i<j} T_ij < min_i e_i.
4. K=5: worst case of a 3-subset majority over P(p,r) for K=5 == K=3 closed form on the subset?
5. Max-min over ALL randomized predictors (BF game, equality constraints) vs randomized self-dual:
   equal when phi symmetric (phi(v)=phi(1-v)); compare for asymmetric phi (true joint law).
"""
import json, itertools
import numpy as np
from scipy.optimize import linprog
from core import *

rng = np.random.default_rng(99)
out = {}

def sharp_vals(vals, err, mu, s):
    A = np.vstack([np.ones_like(vals), vals, vals ** 2])
    res = linprog(-err.astype(float), A_eq=A, b_eq=[1, mu, s], bounds=(0, None), method="highs")
    return (-res.fun) if res.status == 0 else None

def closed_form_K(K, mu, s):
    r = (K + 1) // 2
    m1, m2 = K * mu, K * K * s
    cands = [1.0, ((r + K) * m1 - m2) / (r * K)]
    for j in range(0, r - 1):
        cands.append((m2 - (2 * j + 1) * m1 + j * (j + 1)) / ((r - j) * (r - j - 1)))
    return min(cands)

def closed_form_vals(vals, mu, s):
    """General value set V (sorted, contains 0 and 1), error iff v > 1/2."""
    V = np.unique(vals)
    wstar = V[V > 0.5].min()
    lows = V[V <= 0.5]
    cands = [1.0, ((wstar + 1) * mu - s) / wstar]           # h(w) = w (w* + 1 - w) / w*
    for a, b in zip(lows[:-1], lows[1:]):                      # h(w) = (w-a)(w-b)/((w*-a)(w*-b))
        cands.append((s - (a + b) * mu + a * b) / ((wstar - a) * (wstar - b)))
    return min(cands)

# 1. finite K closed form
maxerr = {}
for K in (3, 5, 7, 11, 21, 101):
    n = np.arange(K + 1)
    vals = n / K
    errs = []
    for _ in range(1500):
        q = rng.dirichlet(np.ones(K + 1) * rng.uniform(0.05, 2))
        mu, s = float(q @ vals), float(q @ vals ** 2)
        lp = sharp_vals(vals, n > K / 2, mu, s)
        errs.append(abs(lp - closed_form_K(K, mu, s)))
    maxerr[K] = float(np.max(errs))
out["finiteK_closed_form_max_err"] = maxerr
print("finite-K closed form vs LP, max abs err:", maxerr)
# continuum polynomial bound validity
w = np.linspace(0, 1, 100001)
out["poly_3w_minus_2w2_minus_indicator_min"] = float(np.min(3 * w - 2 * w ** 2 - (w >= 0.5)))
print("min_w [3w-2w^2 - 1{w>=1/2}] on [0,1] =", out["poly_3w_minus_2w2_minus_indicator_min"])

# 2. weighted
for K in (3, 5):
    Z = patterns(K)
    errs, gaps_cont = [], []
    for _ in range(300):
        pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3), acc_range=(0.55, 0.9))
        wts = rng.exponential(size=K)
        rho = wts / wts.sum()
        W = (1 - Z) @ rho
        if np.any(np.abs(W - 0.5) < 1e-9):
            continue
        mu, s = float(pi @ W), float(pi @ W ** 2)
        lp_full = sharp_vals(W, W > 0.5, mu, s)                    # LP over {0,1}^K (2^K variables)
        cf = closed_form_vals(W, mu, s)
        errs.append(abs(lp_full - cf))
        cont = min((s - mu ** 2) / (s - mu + 0.25) if mu < 0.5 else 1, 3 * mu - 2 * s, 1)
        gaps_cont.append(cont - lp_full)
    out[f"weighted_K{K}_closed_form_max_err"] = float(np.max(errs))
    out[f"weighted_K{K}_continuum_minus_sharp_mean"] = float(np.mean(gaps_cont))
    print(f"weighted K={K}: closed form (value-set) vs 2^K LP max err {np.max(errs):.2e}; continuum bound - sharp: mean {np.mean(gaps_cont):.4f}")

# 3. K=3 identities
K = 3
Z = patterns(K)
e1, e2 = 0.0, 0
G, _ = enumerate_selfdual_ltfs(3, n_random=1000)
for _ in range(2000):
    pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3), acc_range=(0.5, 0.95))
    p = marginals(pi, K); e = 1 - p
    r = pairwise_joint(pi, K)
    T = sum(1 - p[i] - p[j] + r[(i, j)] for (i, j) in r)
    D = 1 - pi[0] - pi[7]
    e1 = max(e1, abs(D - (e.sum() - T)))
    A_p, b_p = set_pairwise(p, r, K)
    lo = [lp_min(g.astype(float), A_p, b_p) for g in G]
    k = int(np.argmax(lo))
    is_maj = G[k].sum() == 4 and all(G[k][v] == (Z[v].sum() >= 2) for v in range(8))
    if abs(T - e.min()) > 1e-9:
        e2 += int(is_maj == (T < e.min()))
    else:
        e2 += 1
out["K3_D_identity_max_err"] = e1
out["K3_tandem_rule_agree_of_2000"] = e2
print(f"K=3: |D - (sum e - sum T)| max {e1:.2e}; 'majority iff sum T < min e' agrees with LP in {e2}/2000")

# 4. subset projection K=5
K = 5
Z = patterns(K)
diffs = []
for _ in range(200):
    pi = mixture_of_products(K, rng, diff_scale=rng.uniform(0, 3), acc_range=(0.55, 0.9))
    p = marginals(pi, K); r = pairwise_joint(pi, K)
    A_p, b_p = set_pairwise(p, r, K)
    S = tuple(sorted(rng.choice(5, 3, replace=False)))
    g = (Z[:, list(S)].sum(1) >= 2).astype(float)
    lo5 = lp_min(g, A_p, b_p)
    e = 1 - p
    T = sum(1 - p[i] - p[j] + r[(i, j)] for i, j in itertools.combinations(S, 2))
    # d_i^S: lone-dissenter prob inside S, from pairwise disagreements within S
    dis = lambda i, j: e[i] + e[j] - 2 * (1 - p[i] - p[j] + r[(min(i, j), max(i, j))])
    dS = [(dis(i, j) + dis(i, k) - dis(j, k)) / 2 for (i, j, k) in ((S[0], S[1], S[2]), (S[1], S[0], S[2]), (S[2], S[0], S[1]))]
    worst_err3 = min(T, min(e[S[a]] + dS[a] for a in range(3)))
    diffs.append(abs((1 - lo5) - worst_err3))
out["K5_subset_majority_vs_K3_closed_form_max_diff"] = float(np.max(diffs))
out["K5_subset_majority_frac_equal"] = float(np.mean(np.array(diffs) < 1e-7))
print(f"K=5: 3-subset majority worst-case (5-var LP) vs K=3 closed form: max diff {np.max(diffs):.2e}, equal in {np.mean(np.array(diffs)<1e-7):.3f}")

# 5. BF vs randomized self-dual
res_sym, res_asym = [], []
for K in (3, 5):
    for _ in range(100 if K == 3 else 60):
        alpha, beta, eta = random_joint(K, rng, diff_scale=rng.uniform(0, 3))
        pi = alpha + beta; p = marginals(pi, K)
        phi = phi_from(alpha, beta, K)
        m = pair_masses(phi, K)
        A_e, b_e = set_even(p, m, K)
        rsd, _ = maxmin_selfdual_randomized(A_e, b_e, K)
        bf_sym, _ = bf_value(0.5 * m, p, K)
        bf_true, _ = bf_value(phi, p, K)
        res_sym.append(bf_sym - rsd)
        res_asym.append(bf_true - rsd)
res_sym, res_asym = np.array(res_sym), np.array(res_asym)
out["BF_sym_minus_randselfdual_maxabs"] = float(np.max(np.abs(res_sym)))
out["BF_truephi_minus_randselfdual_mean"] = float(res_asym.mean())
out["BF_truephi_minus_randselfdual_max"] = float(res_asym.max())
out["BF_truephi_frac_strict_gain"] = float(np.mean(res_asym > 1e-7))
print(f"BF(sym phi) - rand self-dual: max |diff| {np.max(np.abs(res_sym)):.2e}; BF(true asym phi) - rand self-dual: mean {res_asym.mean():.4f}, "
      f"max {res_asym.max():.4f}, strict gain in {np.mean(res_asym>1e-7):.3f}")
json.dump(out, open("results/q3b_results.json", "w"), indent=2)
