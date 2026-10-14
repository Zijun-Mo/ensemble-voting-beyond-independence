"""Figures for the notes.  Run after q1..q4 scripts (reads results/)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from q3_ml_bounds import sharp_K, sharp_cont, c_bound, feasible_K  # noqa (module runs its experiments on import)

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID, BASE, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": BASE, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "legend.frameon": False, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlecolor": INK})

# ---------- Fig 1: sharp (mu,s) bound vs PAC-Bayes oracle bounds ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
mu = 0.2
ss = np.linspace(mu ** 2 + 1e-4, mu - 1e-4, 300)
ax.plot(ss, np.minimum(2 * mu, 1) * np.ones_like(ss), color=C[3], lw=2, label="first-order 2E[L]")
ax.plot(ss, np.minimum(4 * ss, 1), color=C[4], lw=2, label="tandem 4E[L(h,h')]")
ax.plot(ss, np.minimum(c_bound(mu, ss), 1), color=C[1], lw=2, label="C-bound = min_nu Wu et al. Thm 7")
ax.plot(ss, sharp_cont(mu, ss), color=INK, lw=2, ls="--", label="sharp, K -> inf: min(C, 3mu-2s)")
for K, col in ((5, C[0]), (11, C[2])):
    v = np.array([sharp_K(K, mu, s) if feasible_K(K, mu, s) else np.nan for s in ss])
    ax.plot(ss, v, color=col, lw=2, label=f"sharp, K={K} (LP)")
ax.axvline(mu / 2, color=MUTED, lw=1, ls=":")
ax.text(mu / 2 + 0.002, 0.03, "s = mu/2", color=MUTED, fontsize=9)
ax.set_xlabel("s = E_rho2[tandem loss]  (mu = E_rho[L] = 0.2)")
ax.set_ylabel("upper bound on L(MV)")
ax.set_title("Bounds on majority-vote error given (mu, s)")
ax.set_ylim(0, 0.62)
ax.legend(fontsize=8, loc="upper left")
ax = axes[1]
res = json.load(open("results/q3_results.json"))["grid"]
Ks = ["3", "5", "11", "21", "101", "inf"]
ratio = [res[k]["mean_ratio_sharp_over_best_existing"] for k in Ks]
gap = [res[k]["mean_gap_best_existing_minus_sharp"] for k in Ks]
x = np.arange(len(Ks))
ax.bar(x - 0.18, gap, width=0.34, color=C[0], label="mean gap: best existing - sharp")
ax.bar(x + 0.18, ratio, width=0.34, color=C[1], label="mean ratio: sharp / best existing")
for i in range(len(Ks)):
    ax.text(x[i] - 0.18, gap[i] + 0.01, f"{gap[i]:.3f}", ha="center", fontsize=8, color=INK2)
    ax.text(x[i] + 0.18, ratio[i] + 0.01, f"{ratio[i]:.2f}", ha="center", fontsize=8, color=INK2)
ax.set_xticks(x); ax.set_xticklabels(["K=" + k for k in Ks])
ax.set_ylim(0, 1.0)
ax.set_title("Grid over (mu, s): how loose is min(2mu, 4s, C)?")
ax.legend(fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig("figs/fig1_sharp_vs_pacbayes.png", dpi=160)
plt.close(fig)

# ---------- Fig 2: symmetric phase diagram of max-min rule ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 5.4))
q2c = json.load(open("results/q2c_results.json"))
for ax, K in zip(axes, (3, 5)):
    wv = np.load(f"results/q2c_worst_K{K}.npy")
    best = np.full(wv.shape[:2], -1)
    okm = ~np.isnan(wv[..., 0])
    # ties (within 1e-7, e.g. rho = 1 where every rule has worst case p) resolved toward the simplest rule (lowest index)
    mx = np.nanmax(np.where(okm[..., None], wv, -np.inf), axis=2)
    first = np.argmax(wv >= mx[..., None] - 1e-7, axis=2)
    best[okm] = first[okm]
    ps = np.array(q2c[f"K{K}"]["ps"]); rhos = np.array(q2c[f"K{K}"]["rhos"])
    types = q2c[f"K{K}"]["types"]
    used = [k for k in range(len(types)) if np.any(best == k)]
    cmap_cols = ["#f0efec"] + [C[i] for i in range(len(used))]
    remap = np.full(best.shape, 0)
    for j, k in enumerate(used):
        remap[best == k] = j + 1
    ax.pcolormesh(ps, rhos, remap, cmap=ListedColormap(cmap_cols), vmin=-0.5, vmax=len(used) + 0.5, shading="nearest")
    for j, k in enumerate(used):
        ax.plot([], [], "s", color=C[j], ms=9, label=f"weights {types[k]}")
    ax.plot([], [], "s", color="#f0efec", ms=9, label="(p, rho) infeasible")
    if K == 3:
        thr = ((5 * ps - 2) / 3 - ps ** 2) / (ps * (1 - ps))
        ax.plot(ps, thr, color=INK, lw=1.5, ls="--", label="rho*(p) = ((5p-2)/3 - p^2)/(p(1-p))")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("common accuracy p"); ax.set_ylabel("common pairwise correlation rho of Z")
    ax.set_title(f"K={K}: max-min optimal self-dual rule given (p, rho)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2)
    ax.grid(False)
fig.tight_layout()
fig.savefig("figs/fig2_phase_symmetric.png", dpi=160)
plt.close(fig)

# ---------- Fig 3: information hierarchy (K=5) ----------
fig, ax = plt.subplots(figsize=(10, 4.2))
names = ["true error", "sharp | p + full unlabeled law", "sharp | p + pairwise", "sharp | (mu, s)",
         "C-bound", "first-order 2mu", "tandem 4s"]
for j, w in enumerate(("uniform", "logodds")):
    R = np.minimum(np.load(f"results/q3_hierarchy_{w}.npy"), 1)
    m = R.mean(0)
    x = np.arange(len(names)) + (j - 0.5) * 0.38
    ax.bar(x, m, width=0.36, color=C[j], label="uniform weights" if w == "uniform" else "log-odds weights")
    for xi, mi in zip(x, m):
        ax.text(xi, mi + 0.008, f"{mi:.3f}", ha="center", fontsize=8, color=INK2)
ax.set_xticks(np.arange(len(names))); ax.set_xticklabels(names, rotation=15, fontsize=8.5)
ax.set_ylabel("mean over 287 random instances")
ax.set_title("K=5: majority-vote error bounds under increasingly rich information")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("figs/fig3_information_hierarchy.png", dpi=160)
plt.close(fig)

# ---------- Fig 4: K=5 majority intervals, pairwise vs full unlabeled ----------
rows = np.load("results/q1_K5_intervals.npy")
order = np.argsort(rows[:, 6])[::5]
fig, ax = plt.subplots(figsize=(10, 4.2))
xx = np.arange(len(order))
ax.vlines(xx - 0.15, rows[order, 0], rows[order, 1], color=C[1], lw=2, label="sharp interval | p + pairwise agreement")
ax.vlines(xx + 0.15, rows[order, 2], rows[order, 3], color=C[0], lw=2, label="sharp interval | p + full unlabeled law")
ax.plot(xx, rows[order, 6], "o", color=INK, ms=4, label="true majority accuracy")
ax.set_xlabel("instances (every 5th of 300, sorted by true accuracy)")
ax.set_ylabel("accuracy of 5-voter majority")
ax.set_title("What unlabeled data buys beyond pairwise agreement (K=5, class-dependent error laws)")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig("figs/fig4_K5_intervals.png", dpi=160)
plt.close(fig)

# ---------- Fig 5: empirical ----------
S = json.load(open("results/q4_summary.json"))
methods = [("best_single(val)", "best single (val)"), ("majority", "majority"), ("logodds", "log-odds"),
           ("maxmin_pw_plug", "max-min | pairwise (plug-in)"), ("maxmin_full_plug", "max-min | full unlabeled (plug-in)")]
fig, axes = plt.subplots(1, 4, figsize=(13, 3.8), sharey=False)
for ax, (dname, d) in zip(axes, S.items()):
    for j, (key, lab) in enumerate(methods):
        ax.errorbar(j, d[key]["mean"], yerr=d[key]["sd"] / np.sqrt(30), fmt="o", color=C[j], ms=7, capsize=3, label=lab)
    ax.axhline(d["oracle_rule(test)"]["mean"], color=MUTED, ls=":", lw=1)
    ax.set_xticks([]); ax.set_title(dname, fontsize=10)
axes[0].set_ylabel("test accuracy (mean +- s.e., 30 splits)")
axes[0].text(0, S["breast_cancer"]["oracle_rule(test)"]["mean"] + 0.0005, "oracle rule", color=MUTED, fontsize=8)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=8.5)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig("figs/fig5_empirical.png", dpi=160)
plt.close(fig)
print("figures written")
