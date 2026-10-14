"""Exp 14: K=3 caps: exact 9-term dual-vertex formula vs LP (unclipped caps, random + adversarial)."""
import numpy as np, json
from core import states, constraint_system, V
rng = np.random.default_rng(14)
Z3 = states(3); maj = (Z3.sum(1) >= 2).astype(int)
def nine(q, e):
    e12, e13, e23 = e
    return min(1, q.sum() / 2, e12 + e13 + e23, e12 + q[2], e13 + q[1], e23 + q[0], q[0] + q[1], q[0] + q[2], q[1] + q[2])
err = 0; n = 0; crit = 0
for t in range(4000):
    p = rng.uniform(0, 1, 3) if t % 3 == 0 else rng.uniform(0.5, 1, 3); q = 1 - p
    e = rng.uniform(0, 1, 3) * (rng.choice([0.3, 1, 2]))
    rb = np.full((3, 3), np.nan); rb[0, 1] = e[0] - 1 + p[0] + p[1]; rb[0, 2] = e[1] - 1 + p[0] + p[2]; rb[1, 2] = e[2] - 1 + p[1] + p[2]
    v = V(maj, p, rbar=rb)
    if v is None: continue
    n += 1; err = max(err, abs((1 - v) - nine(q, e)))
    if p.min() >= 0.5:
        crit += int((v > p.max() + 1e-12) == (e.sum() < q.min() - 1e-12))
print(dict(n=n, max_abs_err=err, criterion_agree=crit))
json.dump(dict(n=n, max_abs_err=err), open("exp14_results.json", "w"))
