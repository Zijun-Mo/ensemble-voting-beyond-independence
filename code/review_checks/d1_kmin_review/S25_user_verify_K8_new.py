#!/usr/bin/env python3
"""
(Copied verbatim from the user's message of 2026-10-13 23:24.)
Verification of the newly identified K=8 extreme points in ext(Q_{K_8}^uparrow):
  1. Denominator-4 orbit (m = 1 zero-cost edge, c = 7, values in {0, 1/4, 1/2, 3/4})
  2. Two Denominator-5 orbits with m = 0 zero-cost edges (c = 8, values in {1/5, 2/5, 3/5})
All three satisfy rank(A_up8_tight) == 28 and exact B_8(u) = 1.000000.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
import numpy as np
import scipy.optimize as opt

K = 8
edges8 = list(combinations(range(K), 2))
edge_idx8 = {e: i for i, e in enumerate(edges8)}
m8 = len(edges8)

# Build all 3608 odd cycles of K_8 (lengths 3, 5, 7)
cycles8 = []
for length in (3, 5, 7):
    for nodes in combinations(range(K), length):
        first = nodes[0]
        for perm in permutations(nodes[1:]):
            if perm > perm[::-1]:
                continue
            cyc = (first,) + perm
            row = np.zeros(m8)
            for i in range(length):
                u, v = cyc[i], cyc[(i + 1) % length]
                if u > v:
                    u, v = v, u
                row[edge_idx8[(u, v)]] = 1.0
            cycles8.append(row)
A_cyc8 = np.array(cycles8)
A_up8 = np.vstack([A_cyc8, np.eye(m8)])
b_up8 = np.concatenate([np.ones(len(A_cyc8)), np.zeros(m8)])


def compute_exact_B(K, v_frac):
    edges = list(combinations(range(K), 2))
    m_edges = len(edges)
    b_vec = np.array([float(val) for val in v_frac])
    zero_edges = [edges[idx] for idx in range(m_edges) if v_frac[idx] == 0]

    parent = list(range(K))
    parity = [0] * K

    def find(x):
        if parent[x] != x:
            r = find(parent[x])
            parity[x] ^= parity[parent[x]]
            parent[x] = r
        return parent[x]

    for u, w in zero_edges:
        ru, rw = find(u), find(w)
        if ru != rw:
            parent[rw] = ru
            parity[rw] = parity[u] ^ parity[w] ^ 1

    roots = sorted(set(find(x) for x in range(K)))
    comp_map = {r: idx for idx, r in enumerate(roots)}
    verts = [(comp_map[find(x)], parity[x]) for x in range(K)]
    c = len(roots)

    paid_edges = [
        (idx, edges[idx][0], edges[idx][1], b_vec[idx])
        for idx in range(m_edges)
        if v_frac[idx] > 0
    ]
    states = list(product([0, 1], repeat=c - 1))
    m_p = len(paid_edges)
    num_z = len(states)
    total_vars = m_p + num_z
    c_obj = np.zeros(total_vars)
    for i, (_, _, _, cost) in enumerate(paid_edges):
        c_obj[i] = cost
    integrality = np.zeros(total_vars)
    integrality[m_p:] = 1
    A_ub, b_ub = [], []
    for s_idx, bits in enumerate(states):
        spins = (0,) + bits
        z_col = m_p + s_idx
        rowA = np.zeros(total_vars)
        rowD = np.zeros(total_vars)
        for e_idx, (_, u, w, _) in enumerate(paid_edges):
            cu, su = verts[u]
            cw, sw = verts[w]
            in_S_u = su ^ spins[cu]
            in_S_w = sw ^ spins[cw]
            if in_S_u == 1 and in_S_w == 1:
                rowA[e_idx] = -1.0
            elif in_S_u == 0 and in_S_w == 0:
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
        bounds=opt.Bounds(np.zeros(total_vars), np.ones(total_vars) * 10.0),
        constraints=opt.LinearConstraint(np.array(A_ub), -np.inf, np.array(b_ub)),
    )
    return res.fun, c, len(zero_edges)


# 1. Explicit Denominator-4 extreme point in ext(Q_{K_8}^uparrow) (m = 1, c = 7)
u_denom4 = tuple(
    Fraction(s)
    for s in [
        "1/2", "1/4", "1/4", "1/4", "1/4", "1/2", "1/2",
        "3/4", "3/4", "1/4", "1/4", "1/2", "0",   "1/2",
        "1/2", "1/2", "1/4", "1/4", "1/2", "1/2", "1/4",
        "1/4", "1/2", "1/4", "3/4", "1/4", "3/4", "1/2",
    ]
)

# 2. Explicit m = 0, Denominator-5 Orbit I (C_5 + P_3 on K_8)
u_d5_m0_I_dict = {}
for u in range(5):
    for v in range(u + 1, 5):
        d = min((u - v) % 5, (v - u) % 5)
        u_d5_m0_I_dict[(u, v)] = Fraction(1, 5) if d == 1 else Fraction(3, 5)
u_d5_m0_I_dict[(5, 6)] = Fraction(1, 5)
u_d5_m0_I_dict[(5, 7)] = Fraction(1, 5)
u_d5_m0_I_dict[(6, 7)] = Fraction(3, 5)
for u in range(5):
    for w in (5, 6, 7):
        u_d5_m0_I_dict[(u, w)] = Fraction(2, 5)
u_d5_m0_I = tuple(u_d5_m0_I_dict[e] for e in edges8)

# 3. Explicit m = 0, Denominator-5 Orbit II (C_5 with twin vertex 7 suspended by K_2 = (5, 6))
u_d5_m0_II_dict = dict(u_d5_m0_I_dict)
u_d5_m0_II_dict[(5, 6)] = Fraction(1, 5)
for u in range(5):
    u_d5_m0_II_dict[(u, 5)] = u_d5_m0_II_dict[(u, 6)] = Fraction(2, 5)
# Vertex 7 is a twin of vertex 0 on C_5 (with edge (0, 7) having weight 3/5)
u_d5_m0_II_dict[(0, 7)] = Fraction(3, 5)
for u in range(1, 7):
    u_d5_m0_II_dict[(u, 7)] = u_d5_m0_II_dict[(min(0, u), max(0, u))]
u_d5_m0_II = tuple(u_d5_m0_II_dict[e] for e in edges8)

for name, u_vec in [
    ("Denom-4 (m=1, c=7)", u_denom4),
    ("Denom-5 Orbit I (m=0, c=8, C_5 + P_3)", u_d5_m0_I),
    ("Denom-5 Orbit II (m=0, c=8, twin-C_5 + K_2)", u_d5_m0_II),
]:
    flt = np.array([float(x) for x in u_vec])
    slacks = A_up8 @ flt - b_up8
    rk = np.linalg.matrix_rank(A_up8[np.abs(slacks) < 1e-7], tol=1e-6)
    val, c_comp, n0 = compute_exact_B(8, u_vec)
    print(f"{name}: in Q^uparrow={np.min(slacks) >= -1e-7}, tight_rank={rk}/28, #0={n0}, exact B_8(u) = {val:.6f}")

# --- appended by reviewer: save the three vectors for exact cross-checks ---
import pickle, sys
pickle.dump({'denom4': u_denom4, 'd5_I': u_d5_m0_I, 'd5_II': u_d5_m0_II}, open(sys.argv[1], 'wb')) if len(sys.argv) > 1 else None
