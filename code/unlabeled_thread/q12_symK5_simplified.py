"""Q12: symmetric K=5 max-min error simplifies to  min{e, 3T, 1-5e+10T}  (e=1-p, T=pairwise joint error).
Check against the saved LP grid of q2c (worst-case accuracy per rule type, results/q2c_worst_K5.npy)."""
import json, numpy as np
W = np.load("results/q2c_worst_K5.npy"); d = json.load(open("results/q2c_results.json"))
ps, rhos = np.array(d["K5"]["ps"]), np.array(d["K5"]["rhos"])
errs = []; wins = {"dict": 0, "tri": 0, "maj5": 0}
for a, p in enumerate(ps):
    for b, rho in enumerate(rhos):
        if np.isnan(W[b, a]).any():
            continue
        e = 1 - p; T = e * e + rho * p * e
        c = [e, 3 * T, 1 - 5 * e + 10 * T]
        errs.append(abs((1 - min(c)) - np.nanmax(W[b, a])))
        wins[["dict", "tri", "maj5"][int(np.argmin(c))]] += 1
print(dict(n=len(errs), max_abs_err=float(max(errs)), argmin_counts=wins))
json.dump(dict(n=len(errs), max_abs_err=float(max(errs)), argmin_counts=wins), open("results/q12_results.json", "w"))
