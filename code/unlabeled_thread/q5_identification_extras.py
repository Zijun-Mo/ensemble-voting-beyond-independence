"""Q5 (completion session): identification extras.

(a) K=3 width formula: width of the majority-accuracy interval over P(p,r)=E(m,p) equals
      2 * ( min_{z: #errors even} pi(z) + min_{z: #errors odd} pi(z) )
    (one-dimensional direction delta = chi_{123}(s)).  => point identification iff BOTH minima are 0.
(b) Mazzetto et al. (AISTATS 2021) Prop.5 equivalence: worst-case 3-MV error
      sum e - (D12+D13+D23)/2 - 2q,  q = 1/2 max{e_i+e_j-D_ik-D_jk (3 choices), 0}
    equals our min{sum T, min_i (e_i+d_i)}.
(c) Value of information, K=3 and K=5: majority interval width and max-min value under
      p only  |  p + pairwise agreement  |  p + full unlabeled law (case a)  |  + known eta  |
      class-conditional accuracies (a1_i, a0_i) + eta + full unlabeled law.
(d) Sub-committee projection: is the worst case of a 3-subset majority over the K=5 pairwise set equal to the
    K=3 closed form on the triple, also on adversarial families (anti, symmetric negative rho)?
(e) Weighted vs unweighted-subset max-min (Mazzetto's subset selection = max-min over subset majorities):
    how often / by how much does the best of all 81 self-dual LTFs beat the best of the 16 subset majorities?
"""
import json, itertools
import numpy as np
from scipy.optimize import linprog
from core import *
from q2b_regret_and_anticorr import anti_instance

rng = np.random.default_rng(2026)
out = {}


def nerr(K):
    return (1 - patterns(K)).sum(1)


# ---------------------------------------------------------------- (a) K=3 width formula
K = 3
Z = patterns(K); maj = (Z.sum(1) >= 2).astype(float); ne = nerr(K)
errs, zero_w, n = [], 0, 0
for t in range(3000):
    fam = ["independent", "mixture", "bloc"][t % 3]
    pi = family_instance(K, fam, rng) if t % 4 else anti_instance(K, rng)
    p = marginals(pi, K); r = pairwise_joint(pi, K)
    A, b = set_pairwise(p, r, K)
    lo, hi = lp_range(maj, A, b)
    w_formula = 2 * (pi[ne % 2 == 0].min() + pi[ne % 2 == 1].min())
    errs.append(abs((hi - lo) - w_formula)); n += 1
    zero_w += (hi - lo) < 1e-10
out["a_K3_width_formula"] = dict(n=n, max_abs_err=float(max(errs)), frac_zero_width=zero_w / n)
print("(a)", out["a_K3_width_formula"])

# ---------------------------------------------------------------- (b) Mazzetto Prop 5 equivalence
errs = []
for t in range(3000):
    pi = family_instance(3, ["independent", "mixture", "bloc"][t % 3], rng) if t % 4 else anti_instance(3, rng)
    p = marginals(pi, 3); e = 1 - p; r = pairwise_joint(pi, 3)
    D = {k: 1 - (pi @ (Z[:, k[0]] == Z[:, k[1]]).astype(float)) for k in r}
    D12, D13, D23 = D[(0, 1)], D[(0, 2)], D[(1, 2)]
    q = 0.5 * max(e[0] + e[1] - D13 - D23, e[1] + e[2] - D12 - D13, e[0] + e[2] - D12 - D23, 0)
    mazz = e.sum() - (D12 + D13 + D23) / 2 - 2 * q
    T = [(e[i] + e[j] - D[(i, j)]) / 2 for (i, j) in D]
    d = [(D12 + D13 - D23) / 2, (D12 + D23 - D13) / 2, (D13 + D23 - D12) / 2]
    ours = min(sum(T), min(e[i] + d[i] for i in range(3)))
    A, b = set_pairwise(p, r, 3)
    lp = 1 - lp_min(maj, A, b)
    errs.append(max(abs(mazz - ours), abs(ours - lp)))
out["b_mazzetto_prop5_equiv"] = dict(n=3000, max_abs_err=float(max(errs)))
print("(b)", out["b_mazzetto_prop5_equiv"])


# ---------------------------------------------------------------- (c) value of information
def direct_classcond_range(cz, phi, eta, a1, a0, K):
    n = 2 ** K; Zk = patterns(K); fl = flip_index(K)
    rows, rhs = [], []
    for v in range(n):
        row = np.zeros(2 * n); row[v] = 1; row[n + fl[v]] = 1; rows.append(row); rhs.append(phi[v])
    for i in range(K):
        rows.append(np.concatenate([Zk[:, i], np.zeros(n)])); rhs.append(eta * a1[i])
        rows.append(np.concatenate([np.zeros(n), Zk[:, i]])); rhs.append((1 - eta) * a0[i])
    rows.append(np.concatenate([np.ones(n), np.zeros(n)])); rhs.append(eta)
    return lp_range(np.concatenate([cz, cz]), np.array(rows), np.array(rhs))


def p_only_set(p, K):
    Zk = patterns(K)
    return np.vstack([np.ones(2 ** K), Zk.T]), np.concatenate([[1.0], p])


for K in (3, 5):
    Zk = patterns(K)
    G, W = enumerate_selfdual_ltfs(K, n_random=20000)
    Gf = G.astype(float)
    majK = (Zk.sum(1) > K / 2).astype(float)
    recs = []
    ntr = 300 if K == 3 else 150
    for t in range(ntr):
        alpha, beta, eta = random_joint(K, rng, acc_range=(0.6, 0.85), diff_scale=rng.uniform(0, 3))
        pi = alpha + beta; p = marginals(pi, K)
        if np.any(p <= 0.5):
            continue
        phi = phi_from(alpha, beta, K); m = pair_masses(phi, K); r = pairwise_joint(pi, K)
        a1 = (Zk.T @ alpha) / eta; a0 = (Zk.T @ beta) / (1 - eta)
        sets = {"p_only": p_only_set(p, K), "pairwise": set_pairwise(p, r, K), "full": set_even(p, m, K)}
        rec = {"true_maj": float(majK @ pi)}
        for nm, (A, b) in sets.items():
            lo, hi = lp_range(majK, A, b)
            rec[f"w_{nm}"] = hi - lo; rec[f"lo_{nm}"] = lo
            lows = np.array([lp_min(g, A, b) for g in Gf])
            rec[f"mm_{nm}"] = float(lows.max())
            rec[f"mm_true_{nm}"] = float(Gf[int(np.argmax(lows))] @ pi)
        lo, hi = direct_ab_range(majK, phi, p, K, eta=eta); rec["w_full_eta"] = hi - lo
        lo, hi = direct_classcond_range(majK, phi, eta, a1, a0, K); rec["w_classcond"] = hi - lo; rec["lo_classcond"] = lo
        rec["best_single"] = float(p.max())
        recs.append(rec)
    summ = {k: float(np.mean([r[k] for r in recs])) for k in recs[0]}
    summ["n"] = len(recs)
    summ["frac_classcond_zero_width"] = float(np.mean([r["w_classcond"] < 1e-9 for r in recs]))
    summ["frac_pairwise_mm_gt_best_single"] = float(np.mean([r["mm_pairwise"] > r["best_single"] + 1e-9 for r in recs]))
    summ["frac_full_mm_gt_best_single"] = float(np.mean([r["mm_full"] > r["best_single"] + 1e-9 for r in recs]))
    out[f"c_value_of_info_K{K}"] = summ
    print(f"(c) K={K}", {k: round(v, 4) for k, v in summ.items()})

# ---------------------------------------------------------------- (d) sub-committee projection, (e) weighted vs subset MV
K = 5
Zk = patterns(K)
G, W = enumerate_selfdual_ltfs(K, n_random=20000); Gf = G.astype(float)
subset_idx = []
for k in range(len(G)):
    w = canonical_weights(G[k], K)
    if w is not None and set(w) <= {0, 1}:
        subset_idx.append(k)
out["e_n_subset_majorities"] = len(subset_idx)
tri = list(itertools.combinations(range(K), 3))
d_err, e_rec = [], []


def sym_instance(K, rng):
    """exchangeable law with given p, rho via mixture over number of errors (random feasible)."""
    while True:
        q = rng.dirichlet(np.ones(K + 1) * rng.uniform(0.2, 2))
        nu = nerr(K)
        from math import comb
        pi = np.array([q[int(c)] / comb(K, int(c)) for c in nu])
        if marginals(pi, K).min() > 0.5:
            return pi


fams = ["independent", "mixture", "bloc", "anti", "sym"]
for fam in fams:
    gains, diffs, n = [], [], 0
    for t in range(80):
        if fam == "anti":
            pi = anti_instance(K, rng)
        elif fam == "sym":
            pi = sym_instance(K, rng)
        else:
            pi = family_instance(K, fam, rng)
        p = marginals(pi, K)
        if np.any(p <= 0.5):
            continue
        r = pairwise_joint(pi, K); A, b = set_pairwise(p, r, K)
        # (d)
        for T3 in tri:
            g = (Zk[:, list(T3)].sum(1) >= 2).astype(float)
            w5 = lp_min(g, A, b)
            Z3 = patterns(3); pi3 = np.zeros(8)
            idx3 = (Zk[:, list(T3)] * (2 ** np.arange(3))).sum(1).astype(int)
            np.add.at(pi3, idx3, pi)
            A3, b3 = set_pairwise(marginals(pi3, 3), pairwise_joint(pi3, 3), 3)
            w3 = lp_min((Z3.sum(1) >= 2).astype(float), A3, b3)
            d_err.append(w5 - w3)
        # (e) pairwise and full
        lows = np.array([lp_min(g, A, b) for g in Gf])
        m = pair_masses(pi, K)  # phi -> m: under case (a) m = pi + pi~
        Af, bf = set_even(p, m, K)
        lows_f = np.array([lp_min(g, Af, bf) for g in Gf])
        gains.append((lows.max() - lows[subset_idx].max(), lows_f.max() - lows_f[subset_idx].max()))
        n += 1
    g = np.array(gains)
    e_rec.append(dict(fam=fam, n=n, pw_frac_strict_gain=float(np.mean(g[:, 0] > 1e-7)), pw_mean_gain=float(g[:, 0].mean()),
                      pw_max_gain=float(g[:, 0].max()), full_frac_strict_gain=float(np.mean(g[:, 1] > 1e-7)),
                      full_mean_gain=float(g[:, 1].mean()), full_max_gain=float(g[:, 1].max())))
    print("(e)", e_rec[-1])
d_err = np.array(d_err)
out["d_subcommittee_projection"] = dict(n=len(d_err), max_abs_diff=float(np.abs(d_err).max()),
                                        frac_K5_worst_gt_K3_worst=float(np.mean(d_err > 1e-7)), min_diff=float(d_err.min()))
print("(d)", out["d_subcommittee_projection"])
out["e_weighted_vs_subset_mv"] = e_rec
json.dump(out, open("results/q5_results.json", "w"), indent=2)
