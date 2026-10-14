#!/usr/bin/env python3
"""K=5 closed-form finite-sample certificate C_248 (Theorem 5.4 of the revised task-1 report).

Input: rational confidence endpoints q_minus, q_plus (length 5) and pairwise upper bounds b (length 10, edge order
itertools.combinations(range(5), 2)), already outward-rounded to rationals.
Output: the certificate value C_248 (exact Fraction), the chosen threshold rule (integer weights w; the rule votes
"error" iff 2*sum_i w_i x_i > sum_i w_i), and the certificate (a, beta, gamma) that proves it.

Every one of the 248 certificates is re-checked here in exact arithmetic (Q(x) >= F_w(x) for all 32 states,
gamma >= 0), so the returned value is a rigorous upper bound on the worst-case error of the chosen rule over the
whole empirical ambiguity set; no MILP and no optimality certificate are needed.
Usage: python3 task1_C248_closed_form.py   (runs a demo); or import and call c248(q_minus, q_plus, b)."""
import json, os, itertools
from fractions import Fraction as F
E = list(itertools.combinations(range(5), 2))
_CERTS = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'U248_certificates.json')))
def _check(c):
    a = F(c['a']); be = [F(x) for x in c['beta']]; ga = [F(x) for x in c['gamma']]; w = c['weights']
    assert all(g >= 0 for g in ga) and sum(w) % 2 == 1
    for m in range(32):
        x = [(m >> i) & 1 for i in range(5)]
        Q = a + sum(b_ * xi for b_, xi in zip(be, x)) + sum(g * x[i] * x[j] for g, (i, j) in zip(ga, E))
        assert Q >= int(2 * sum(wi * xi for wi, xi in zip(w, x)) > sum(w))
    return a, be, ga
CERTS = [(c['label'], c['weights']) + _check(c) for c in _CERTS]
assert len(CERTS) == 248
def c248(q_minus, q_plus, b):
    ql = [F(x) for x in q_minus]; qu = [F(x) for x in q_plus]; bb = [F(x) for x in b]
    assert all(0 <= l <= u for l, u in zip(ql, qu)) and all(x >= 0 for x in bb)
    best = None
    for label, w, a, be, ga in CERTS:
        val = a + sum(bi * qu[i] if bi > 0 else bi * ql[i] for i, bi in enumerate(be)) + sum(g * x for g, x in zip(ga, bb))
        if best is None or val < best[0]:
            best = (val, label, w, (a, be, ga))
    val, label, w, cert = best
    return {'C248': min(val, F(1)), 'rule_weights': w, 'certificate_label': label,
            'a': cert[0], 'beta': cert[1], 'gamma': cert[2]}
if __name__ == '__main__':
    r = c248([F(3, 10)] * 5, [F(3, 10)] * 5, [F(1, 20)] * 10)
    print('symmetric demo q=0.3, b=0.05: C248 =', r['C248'], 'rule weights', r['rule_weights'], '(', r['certificate_label'], ')')
    r = c248([F(1, 5) - F(1, 50)] * 5, [F(1, 5) + F(1, 50)] * 5, [F(1, 25) + F(1, 50)] * 10)
    print('demo q=0.2+-0.02, b=0.06:   C248 =', r['C248'], 'rule weights', r['rule_weights'], '(', r['certificate_label'], ')')
