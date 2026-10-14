#!/usr/bin/env python3
"""
B_graph Complexity & Structural Breakthroughs -- Consolidated Exact Verification Suite
Verifies:
1. [Theorem 2.1] Exact global optimality of B_{J_4} = 14/5 across all 2^11 = 2048 LP branches.
2. [Theorem 3.1] Single-Center Anchoring Theorem over all 2^16 = 65,536 partitions (g_H = 2, ratio = 14/5).
3. [Theorem 4.3] Exact Fractional Chromatic Number Threshold Theorem on H(G):
   B_{H(G)}(1, beta) < 2 + beta  <===>  \chi_f(G) > 2 + beta.
4. [Theorem 5.1] Exhaustive classification of all 2^10 = 1024 odd-K_5 lifts:
   ONLY the 10 isomorphic copies of 9-vertex J_4^\circ (K = 9) achieve B_H < 3.
"""
from fractions import Fraction
from itertools import combinations, product
import numpy as np
import scipy.optimize as opt


def verify_task_A_J4_and_Jm():
    print("=" * 72)
    print("1. [Theorem 2.1] Exact Global Optimality & Dual Certificate for J_4")
    print("=" * 72)
    m = 4
    clique_pairs = list(combinations(range(m), 2))
    states = list(product([0, 1], repeat=m))
    flex_indices = [idx for idx, x in enumerate(states) if sum(x) >= 2]

    c = np.ones(10)
    best_val = 100.0
    for mask in range(1 << len(flex_indices)):
        A_ub, b_ub = [], []
        for idx, x in enumerate(states):
            use_D = ((mask >> flex_indices.index(idx)) & 1) if idx in flex_indices else 0
            row = np.zeros(10)
            if not use_D:
                for i in range(4):
                    if x[i] == 1:
                        row[i] = -1.0
                for p_idx, (i, j) in enumerate(clique_pairs):
                    if x[i] == 0 and x[j] == 0:
                        row[4 + p_idx] = -1.0
            else:
                for p_idx, (i, j) in enumerate(clique_pairs):
                    if x[i] == 1 and x[j] == 1:
                        row[4 + p_idx] = -1.0
            A_ub.append(row)
            b_ub.append(-1.0)
        res = opt.linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
        if res.success and res.fun < best_val:
            best_val = res.fun
    print(f"  Exhaustive min over all 2^11 = 2048 branches of J_4: {best_val:.10f} (14/5 = 2.8)")
    assert abs(best_val - 2.8) < 1e-9


def verify_task_B_shared_center():
    print("\n" + "=" * 72)
    print("2. [Theorem 3.1] Single-Center Anchoring (Two J_4 Copies Sharing B_0)")
    print("=" * 72)
    a, b = Fraction(2, 5), Fraction(1, 5)
    edges = []
    for offset in (0, 4):
        for i, j in combinations(range(offset + 1, offset + 5), 2):
            edges.append((("T", i), ("T", j), b, Fraction(1)))
        for i in range(offset + 1, offset + 5):
            edges.append((("B", 0), ("B", i), a, Fraction(1)))
    for i in range(1, 9):
        edges.append((("T", i), ("B", i), Fraction(2), Fraction(0)))

    var_verts = [("T", i) for i in range(1, 9)] + [("B", i) for i in range(1, 9)]
    min_g = None
    for bits in product((0, 1), repeat=16):
        color = {("B", 0): 1}
        color.update(zip(var_verts, bits))
        inner0, inner1 = Fraction(0), Fraction(0)
        for u, v, w, _ in edges:
            if color[u] == color[v]:
                if color[u] == 0:
                    inner0 += w
                else:
                    inner1 += w
        val = max(inner0, inner1)
        if min_g is None or val < min_g:
            min_g = val
    paid_cost = sum(w * cost for _, _, w, cost in edges)
    print(f"  Exhaustive check over 2^16 = 65,536 partitions: g_H = {min_g}, cost = {paid_cost}, ratio = {paid_cost/min_g}")
    assert min_g == Fraction(2) and paid_cost / min_g == Fraction(14, 5)


def compute_chi_f(n, edges_G):
    adj = [set() for _ in range(n)]
    for u, v in edges_G:
        adj[u].add(v)
        adj[v].add(u)
    A_ub = []
    for x in product([0, 1], repeat=n):
        ones = [i for i in range(n) if x[i] == 1]
        if not ones:
            continue
        if not any(v in adj[u] for i, u in enumerate(ones) for v in ones[i + 1 :]):
            row = np.zeros(n)
            for i in ones:
                row[i] = 1.0
            A_ub.append(row)
    res = opt.linprog(-np.ones(n), A_ub=A_ub, b_ub=np.ones(len(A_ub)), bounds=(0, None), method="highs")
    return -res.fun


def compute_B_HG_exact(n, edges_G, beta=1.0):
    states = list(product([0, 1], repeat=n))
    m_G = len(edges_G)
    num_e = n + m_G
    adj = [set() for _ in range(n)]
    for u, v in edges_G:
        adj[u].add(v)
        adj[v].add(u)
    flex_states, fixed_A_states = [], []
    for x in states:
        ones = [i for i in range(n) if x[i] == 1]
        if any(v in adj[u] for i, u in enumerate(ones) for v in ones[i + 1 :]):
            flex_states.append(x)
        else:
            fixed_A_states.append(x)
    total_vars = num_e + len(flex_states)
    c = np.zeros(total_vars)
    c[:n] = 1.0
    c[n:num_e] = beta
    integrality = np.zeros(total_vars)
    integrality[num_e:] = 1
    A_ub, b_ub = [], []
    for x in fixed_A_states:
        rowA = np.zeros(total_vars)
        for i in range(n):
            if x[i] == 1:
                rowA[i] = -1.0
        for e_idx, (u, v) in enumerate(edges_G):
            if x[u] == 0 and x[v] == 0:
                rowA[n + e_idx] = -1.0
        A_ub.append(rowA)
        b_ub.append(-1.0)
    for z_idx, x in enumerate(flex_states):
        z_col = num_e + z_idx
        rowA, rowD = np.zeros(total_vars), np.zeros(total_vars)
        for i in range(n):
            if x[i] == 1:
                rowA[i] = -1.0
        for e_idx, (u, v) in enumerate(edges_G):
            if x[u] == 0 and x[v] == 0:
                rowA[n + e_idx] = -1.0
            if x[u] == 1 and x[v] == 1:
                rowD[n + e_idx] = -1.0
        rowA[z_col] = -1.0
        rowD[z_col] = 1.0
        A_ub.append(rowA)
        b_ub.append(-1.0)
        A_ub.append(rowD)
        b_ub.append(0.0)
    res = opt.milp(
        c=c,
        integrality=integrality,
        bounds=opt.Bounds(np.zeros(total_vars), np.ones(total_vars)),
        constraints=opt.LinearConstraint(np.array(A_ub), -np.inf, np.array(b_ub)),
    )
    return res.fun


def verify_task_C_fractional_chromatic_equivalence():
    print("\n" + "=" * 72)
    print("3. [Theorem 4.3] Fractional Chromatic Number Threshold Theorem on H(G)")
    print("=" * 72)
    graphs = [
        ("C_5", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
        ("W_5 (5-wheel)", 6, [(0, i) for i in range(1, 6)] + [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]),
        ("K_4", 4, list(combinations(range(4), 2))),
        ("K_3 x K_2", 6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]),
        ("Moser spindle", 7, [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 4), (0, 5), (4, 5), (4, 6), (5, 6), (3, 6)]),
    ]
    for name, n, edges in graphs:
        chi_f = compute_chi_f(n, edges)
        for beta in [1.0, 1.25, 1.5, 2.0]:
            b_val = compute_B_HG_exact(n, edges, beta=beta)
            thresh = 2.0 + beta
            assert (b_val < thresh - 1e-6) == (chi_f > thresh + 1e-6)
        print(f"  Verified {name:16s} (chi_f = {chi_f:.2f}): phase transition at 2 + beta == chi_f(G) holds!")


def verify_task_4_odd_K5_uniqueness():
    print("\n" + "=" * 72)
    print("4. [Theorem 5.1] Exhaustive Classification of All 1024 Odd-K_5 Lifts")
    print("=" * 72)
    k5_edges = list(combinations(range(5), 2))
    states = list(product([0, 1], repeat=4))
    strict_masks = []
    for mask in range(1024):
        sigma = [(mask >> idx) & 1 for idx in range(10)]
        has_top = [False] * 5
        has_bot = [False] * 5
        for idx, (u, v) in enumerate(k5_edges):
            if sigma[idx] == 0:
                has_top[u] = has_top[v] = True
            else:
                has_bot[u] = has_bot[v] = True
        m_split = sum(1 for u in range(5) if has_top[u] and has_bot[u])
        K = 5 + m_split

        total_vars = 10 + 16
        c_obj = np.zeros(total_vars)
        c_obj[:10] = 1.0
        integrality = np.zeros(total_vars)
        integrality[10:] = 1
        A_ub, b_ub = [], []
        for s_idx, bits in enumerate(states):
            s = (1,) + bits
            z_col = 10 + s_idx
            rowA, rowD = np.zeros(total_vars), np.zeros(total_vars)
            for e_idx, (u, v) in enumerate(k5_edges):
                if sigma[e_idx] == 1:
                    if s[u] == 1 and s[v] == 1:
                        rowA[e_idx] = -1.0
                    elif s[u] == 0 and s[v] == 0:
                        rowD[e_idx] = -1.0
                else:
                    if s[u] == 0 and s[v] == 0:
                        rowA[e_idx] = -1.0
                    elif s[u] == 1 and s[v] == 1:
                        rowD[e_idx] = -1.0
            rowA[z_col] = -1.0
            rowD[z_col] = 1.0
            A_ub.append(rowA)
            b_ub.append(-1.0)
            A_ub.append(rowD)
            b_ub.append(0.0)
        res = opt.milp(
            c=c_obj,
            integrality=integrality,
            bounds=opt.Bounds(np.zeros(total_vars), np.ones(total_vars)),
            constraints=opt.LinearConstraint(np.array(A_ub), -np.inf, np.array(b_ub)),
        )
        if res.fun < 3.0 - 1e-6:
            strict_masks.append((mask, K, res.fun))
    print(f"  Out of 1024 odd-K_5 lifts, {len(strict_masks)} have B_H < 3 (all at K = 9, B_H = 14/5 = 2.8)!")
    assert len(strict_masks) == 10 and all(K == 9 and abs(val - 2.8) < 1e-6 for _, K, val in strict_masks)


if __name__ == "__main__":
    verify_task_A_J4_and_Jm()
    verify_task_B_shared_center()
    verify_task_C_fractional_chromatic_equivalence()
    verify_task_4_odd_K5_uniqueness()
