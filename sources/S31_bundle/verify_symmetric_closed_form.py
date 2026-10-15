#!/usr/bin/env python3
"""
Exact Rational & Numerical Reproducibility Suite for:
"Closed-Form Minimax Ensemble Voting under Symmetric Parameters for Arbitrary K"

Verifies:
1. Pointwise validity (Q(x) >= F(x), gamma >= 0) of all three certificate families:
   - Odd-subset Boros-Prekopa certificates B_{2s+1, t}
   - Even-K symmetric layer certificates B^{even}_{2m, t} (dual feasibility for both Maj_{2m-1} and w=(2,1,...,1))
   - Even-K broken-tie anchor certificates B^{anch}_{2m, t} (valid for w=(2,1,...,1))
2. Exact symbolic/rational verification of the three Algebraic Dominance Identities (Theorems 3 & 4).
3. Exact rational witness points (fractions.Fraction) proving irredundancy of every active term for K = 3..10.
4. Exact lower-envelope cell vertex enumeration and symmetry-broken MILP audit for K = 3..10.
"""

import time
from fractions import Fraction as F
from itertools import combinations
import numpy as np
from scipy.optimize import linprog, milp, Bounds, LinearConstraint

def term_dict():
    return ("dict", F(0), F(1), F(0))

def term_tri():
    return ("B_3,0", F(0), F(0), F(3))

def term_odd_bp(s, t):
    u = 2 * s + 1
    denom = (s + 1 - t) * (s - t)
    return (f"B_{{{u},{t}}}", F(t * (t + 1), denom), F(-2 * t * u, denom), F(2 * s * u, denom))

def term_even_layer(m, t):
    u = 2 * m
    denom = (m - t) * (m - 1 - t)
    return (f"B^{{even}}_{{{u},{t}}}", F(t * (t + 1), denom), F(-2 * t * u, denom), F(u * (u - 1), denom))

def term_even_anchor(m, t):
    u = 2 * m
    denom = (m + 1 - t) * (m - t)
    return (f"B^{{anch}}_{{{u},{t}}}", F(t * (t + 1), denom), F(2 - 2 * t * (u + 1), denom), F(2 * (u - 1) * (m + 1), denom))

def eval_term(term, q, b):
    _, a, B, C = term
    return a + B * q + C * b

def b_min_sym(K, q):
    return max(F(2 * r * K * q - r * (r + 1), K * (K - 1)) for r in range(K // 2 + 1))

def get_candidate_terms(K):
    terms = [term_dict(), term_tri()]
    if K % 2 == 1:
        m = (K - 1) // 2
        for s in range(2, m + 1):
            for t in range(1, s):
                terms.append(term_odd_bp(s, t))
    else:
        m = K // 2
        for s in range(2, m):
            for t in range(1, s):
                terms.append(term_odd_bp(s, t))
        for t in range(1, m - 1):
            terms.append(term_even_layer(m, t))
            terms.append(term_even_anchor(m, t))
    return terms

def verify_pointwise_validity(max_K=12):
    print("=== [Check 1] Pointwise Validity of Certificate Polynomials ===")
    for K in range(3, max_K + 1):
        if K % 2 == 1:
            s = (K - 1) // 2
            for t in range(0, s):
                denom = F((s + 1 - t) * (s - t))
                for n in range(K + 1):
                    q_val = F((n - t) * (n - t - 1)) / denom
                    f_val = 1 if n >= s + 1 else 0
                    assert q_val >= f_val
        else:
            m = K // 2
            for t in range(1, m - 1):
                denom1 = F((m - t) * (m - 1 - t))
                denom2 = F((m + 1 - t) * (m - t))
                for x0 in (0, 1):
                    for n_rest in range(2 * m):
                        n_tot = x0 + n_rest
                        w_tot = 2 * x0 + n_rest
                        f_w = 1 if w_tot >= m + 1 else 0
                        f_maj = 1 if n_rest >= m else 0
                        q_even = F((n_tot - t) * (n_tot - t - 1)) / denom1
                        q_anch = F((w_tot - t) * (w_tot - t - 1)) / denom2
                        assert q_even >= f_w and q_even >= f_maj
                        assert q_anch >= f_w
    print(f"  PASS: All certificate polynomials verified pointwise for K = 3..{max_K}.")

def verify_algebraic_identities():
    print("\n=== [Check 2] Exact Algebraic Dominance Identities (Theorems 3 & 4) ===")
    for s in range(3, 15):
        for t in range(1, s - 1):
            for q in [F(1, 7), F(2, 9), F(3, 8), F(5, 11)]:
                for b in [F(1, 30), F(1, 15), F(1, 9)]:
                    val_anch = eval_term(term_even_anchor(s, t), q, b)
                    val_odd_t = eval_term(term_odd_bp(s, t), q, b)
                    val_odd_tp1 = eval_term(term_odd_bp(s, t + 1), q, b)
                    assert val_anch - val_odd_t == F(2, (s + 1 - t) * (s - t)) * (q - b)

                    val_even = eval_term(term_even_layer(s, t), q, b)
                    c1 = F(1) - F(t, 2 * s + 1) + F(2 * s, (2 * s + 1) * (s - 1 - t))
                    c2 = F(t, 2 * s + 1)
                    assert val_even == c1 * val_odd_t + c2 * val_odd_tp1

                    val_prev_odd = eval_term(term_odd_bp(s - 1, t), q, b)
                    L1 = F((2 * s - 2 * t - 3) * (s - t + 1), (2 * s + 1) * (s - t - 1))
                    L2 = F(2 * t, 2 * s + 1)
                    L3 = F(2, 3 * (s - t) * (s - t - 1))
                    assert L1 > 0 and L2 > 0 and L3 > 0
                    assert L1 + L2 + L3 - 1 == F(2 * (3 * t + 1 - s), 3 * (2 * s + 1) * (s - t) * (s - t - 1))
                    assert val_prev_odd == L1 * val_odd_t + L2 * val_odd_tp1 + L3 * (3 * b)
    print("  PASS: Theorems 3(a), 3(b), and 4 verified in exact rational arithmetic for s = 3..14.")

def find_irredundant_terms_with_witnesses(K):
    terms = get_candidate_terms(K)
    uniq = {}
    for name, a, B, C in terms:
        uniq.setdefault((a, B, C), []).append(name)
    u_list = [(names[0], a, B, C) for (a, B, C), names in uniq.items()]
    active = []
    for idx, (name, a0, B0, C0) in enumerate(u_list):
        c = [0.0, 0.0, -1.0]
        A_ub = [[-1.0, 1.0, 0.0]]
        b_ub = [0.0]
        for r in range(K // 2 + 1):
            A_ub.append([2.0 * r / (K - 1), -1.0, 0.0])
            b_ub.append(r * (r + 1) / float(K * (K - 1)))
        for j, (n2, a, B, C) in enumerate(u_list):
            if j == idx:
                continue
            A_ub.append([float(B0 - B), float(C0 - C), 1.0])
            b_ub.append(float(a - a0))
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(0, 0.5), (0, 0.5), (None, 1)], method="highs")
        if res.success and -res.fun > 1e-9:
            q_fl, b_fl = res.x[0], res.x[1]
            witness = None
            for denom_lim in (20, 50, 100, 500, 2000):
                qw = F(q_fl).limit_denominator(denom_lim)
                bw = F(b_fl).limit_denominator(denom_lim)
                if not (F(0) <= qw <= F(1, 2) and b_min_sym(K, qw) <= bw <= qw):
                    continue
                val0 = a0 + B0 * qw + C0 * bw
                other_min = min(a + B * qw + C * bw for j, (_, a, B, C) in enumerate(u_list) if j != idx)
                if val0 < other_min:
                    witness = (qw, bw, val0, other_min - val0)
                    break
            assert witness is not None
            active.append((name, a0, B0, C0, witness))
    return active

def exact_milp_sym_fast(K, q, b, time_limit=30):
    edges = list(combinations(range(K), 2))
    m = len(edges)
    reps = [mask for mask in range(1 << (K - 1))]
    n_reps = len(reps)
    n_cont = 1 + K + m
    obj = np.r_[1.0, np.full(K, float(q)), np.full(m, float(b)), np.zeros(n_reps)]
    lower = np.r_[np.full(1 + K, -np.inf), np.zeros(m + n_reps)]
    upper = np.r_[np.full(n_cont, np.inf), np.ones(n_reps)]
    integrality = np.r_[np.zeros(n_cont), np.ones(n_reps)]
    rows, rhs = [], []
    full_mask = (1 << K) - 1
    for idx, mask in enumerate(reps):
        x = np.array([(mask >> i) & 1 for i in range(K)], dtype=float)
        r1 = np.r_[1.0, x, [x[i] * x[j] for i, j in edges], np.zeros(n_reps)]
        r1[n_cont + idx] = -1.0
        rows.append(r1); rhs.append(0.0)
        xc = 1.0 - x
        r2 = np.r_[1.0, xc, [xc[i] * xc[j] for i, j in edges], np.zeros(n_reps)]
        r2[n_cont + idx] = 1.0
        rows.append(r2); rhs.append(1.0)
    for i in range(K - 1):
        row = np.zeros(n_cont + n_reps)
        row[1 + i] = 1.0
        row[1 + i + 1] = -1.0
        for k, (u, v) in enumerate(edges):
            if u == i or v == i: row[1 + K + k] += 0.5
            if u == i + 1 or v == i + 1: row[1 + K + k] -= 0.5
        rows.append(row); rhs.append(0.0)
    seen_pairs = set()
    for mask in range(1 << K):
        for i in range(K - 1):
            if ((mask >> i) & 1) == 0 and ((mask >> (i + 1)) & 1) == 1:
                mask2 = mask ^ (1 << i) ^ (1 << (i + 1))
                c1 = mask if (mask >> (K - 1)) == 0 else (full_mask ^ mask)
                c2 = mask2 if (mask2 >> (K - 1)) == 0 else (full_mask ^ mask2)
                key = (min(c1, c2), max(c1, c2), (mask2 >> (K - 1)) != (mask >> (K - 1)))
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)
                row = np.zeros(n_cont + n_reps)
                rhs_val = 0.0
                if (mask2 >> (K - 1)) == 0:
                    row[n_cont + mask2] += 1.0
                else:
                    row[n_cont + (full_mask ^ mask2)] -= 1.0
                    rhs_val -= 1.0
                if (mask >> (K - 1)) == 0:
                    row[n_cont + mask] -= 1.0
                else:
                    row[n_cont + (full_mask ^ mask)] += 1.0
                    rhs_val += 1.0
                rows.append(row); rhs.append(rhs_val)
    ans = milp(obj, integrality=integrality, bounds=Bounds(lower, upper),
               constraints=LinearConstraint(np.array(rows), np.array(rhs), np.inf),
               options={"mip_rel_gap": 1e-12, "time_limit": time_limit})
    return float(ans.fun), ans.mip_gap

def verify_vertices(max_K=10):
    print("\n=== [Check 3 & 4] Irredundant Terms, Rational Witnesses & Vertex Audits ===")
    for K in range(3, max_K + 1):
        active = find_irredundant_terms_with_witnesses(K)
        print(f"\n--- K = {K} ({len(active)} irredundant terms) ---")
        for name, a, B, C, (qw, bw, val0, margin) in active:
            print(f"  Term {name:18s}: L(q,b) = {str(a):>5s} + ({str(B):>6s})*q + ({str(C):>6s})*b | witness (q={qw}, b={bw}) -> val={val0}, margin={margin}")
        terms = [(name, a, B, C) for name, a, B, C, _ in active]
        lines = [(F(1), F(0), F(0)), (F(1), F(0), F(1, 2)), (F(1), F(-1), F(0))]
        for r in range(K // 2 + 1):
            lines.append((F(2 * r * K), F(-K * (K - 1)), F(r * (r + 1))))
        for (_, a1, B1, C1), (_, a2, B2, C2) in combinations(terms, 2):
            if B1 != B2 or C1 != C2:
                lines.append((B1 - B2, C1 - C2, a2 - a1))
        pts = set()
        for (A1, B1, D1), (A2, B2, D2) in combinations(lines, 2):
            det = A1 * B2 - A2 * B1
            if det == 0:
                continue
            q = (D1 * B2 - D2 * B1) / det
            b = (A1 * D2 - A2 * D1) / det
            if 0 <= q <= F(1, 2) and b <= q:
                bmin = b_min_sym(K, q)
                if b >= bmin:
                    u_val = min(a + B * q + C * b for _, a, B, C in terms)
                    n_tight_terms = sum(1 for _, a, B, C in terms if a + B * q + C * b == u_val)
                    n_tight_bnds = (q == 0) + (q == F(1, 2)) + (b == q) + sum(
                        1 for r in range(K // 2 + 1) if F(2 * r * K * q - r * (r + 1), K * (K - 1)) == b
                    )
                    if n_tight_terms >= 3 or (n_tight_terms >= 2 and n_tight_bnds >= 1):
                        pts.add((q, b))
        nontriv = [
            (q, b) for q, b in sorted(pts)
            if min(a + B * q + C * b for _, a, B, C in terms) > 0 and q < F(1, 2)
        ]
        print(f"  Nontrivial interior/boundary-tie vertices for K={K}: {len(nontriv)}")
        for q, b in nontriv:
            u_val = min(a + B * q + C * b for _, a, B, C in terms)
            tight = [name for name, a, B, C in terms if a + B * q + C * b == u_val]
            v_milp, gap = exact_milp_sym_fast(K, q, b)
            diff = abs(v_milp - float(u_val))
            assert diff < 1e-9
            print(f"    Vertex (q={str(q):>6s}, b={str(b):>7s}): U_det = {str(u_val):>6s} ({float(u_val):.6f}), MILP = {v_milp:.6f}, diff = {diff:.1e}, tight = {tight}")

if __name__ == "__main__":
    t_start = time.time()
    verify_pointwise_validity(max_K=12)
    verify_algebraic_identities()
    verify_vertices(max_K=10)
    print(f"\nALL REPRODUCIBILITY CHECKS PASSED in {time.time() - t_start:.2f}s.")
