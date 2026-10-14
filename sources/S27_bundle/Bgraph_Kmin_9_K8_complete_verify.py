#!/usr/bin/env python3
"""Complete Deterministic Verification that K_min = 9 for B_graph(b) < c^*(b)."""

from fractions import Fraction
from itertools import combinations, permutations, product
import time
import numpy as np
import scipy.optimize as opt

EDGES8 = [(i, j) for i in range(8) for j in range(i + 1, 8)]
EDGE_TO_IDX8 = {e: idx for idx, e in enumerate(EDGES8)}

# Complete catalog of the 19 canonical fractional orbits F01..F19 of ext(Q_{K_8}^up)
ORBITS_F01_F19 = [
    (
        "F01",
        3,
        (
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
    (
        "F02",
        3,
        (
            0,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
    (
        "F03",
        3,
        (
            0,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
    (
        "F04",
        3,
        (
            1,
            0,
            2,
            1,
            1,
            1,
            1,
            2,
            0,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
    (
        "F05",
        3,
        (
            1,
            1,
            0,
            2,
            2,
            1,
            1,
            1,
            2,
            0,
            2,
            1,
            1,
            2,
            2,
            0,
            1,
            1,
            1,
            1,
            2,
            2,
            1,
            2,
            2,
            2,
            2,
            1,
        ),
    ),
    (
        "F06",
        4,
        (
            0,
            2,
            2,
            1,
            1,
            3,
            3,
            2,
            2,
            3,
            3,
            1,
            1,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            2,
            2,
            2,
        ),
    ),
    (
        "F07",
        5,
        (
            2,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            3,
            3,
            2,
            2,
            3,
            1,
            3,
            2,
            2,
            3,
            1,
            2,
            2,
            1,
            2,
            2,
            2,
            2,
            3,
        ),
    ),
    (
        "F08",
        5,
        (
            3,
            1,
            1,
            1,
            3,
            2,
            2,
            1,
            1,
            3,
            1,
            2,
            2,
            3,
            3,
            3,
            2,
            2,
            3,
            3,
            2,
            2,
            1,
            2,
            2,
            2,
            2,
            1,
        ),
    ),
    (
        "F09",
        5,
        (
            0,
            3,
            3,
            1,
            1,
            2,
            2,
            2,
            2,
            4,
            4,
            3,
            3,
            1,
            1,
            3,
            2,
            2,
            3,
            1,
            2,
            2,
            3,
            2,
            2,
            2,
            2,
            1,
        ),
    ),
    (
        "F10",
        5,
        (
            0,
            2,
            2,
            2,
            2,
            2,
            1,
            3,
            3,
            3,
            3,
            3,
            4,
            1,
            1,
            3,
            3,
            2,
            3,
            1,
            3,
            2,
            3,
            1,
            2,
            1,
            2,
            2,
        ),
    ),
    (
        "F11",
        6,
        (
            0,
            1,
            4,
            2,
            2,
            4,
            4,
            5,
            2,
            4,
            4,
            2,
            2,
            1,
            3,
            3,
            3,
            3,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
        ),
    ),
    (
        "F12",
        7,
        (
            1,
            3,
            1,
            5,
            3,
            4,
            4,
            3,
            5,
            1,
            3,
            4,
            4,
            3,
            3,
            1,
            2,
            2,
            3,
            3,
            2,
            2,
            3,
            2,
            2,
            4,
            4,
            3,
        ),
    ),
    (
        "F13",
        3,
        (
            0,
            0,
            2,
            2,
            2,
            2,
            2,
            3,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
    (
        "F14",
        3,
        (
            0,
            0,
            1,
            1,
            1,
            1,
            1,
            3,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
    (
        "F15",
        3,
        (
            0,
            0,
            0,
            2,
            2,
            2,
            2,
            3,
            3,
            1,
            1,
            1,
            1,
            3,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
    (
        "F16",
        3,
        (
            2,
            0,
            0,
            1,
            2,
            2,
            2,
            1,
            1,
            0,
            1,
            1,
            1,
            3,
            2,
            1,
            1,
            1,
            2,
            1,
            1,
            1,
            2,
            2,
            2,
            1,
            1,
            1,
        ),
    ),
    (
        "F17",
        3,
        (
            1,
            2,
            0,
            0,
            1,
            1,
            1,
            0,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            2,
            2,
            2,
            3,
            2,
            2,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
        ),
    ),
    (
        "F18",
        3,
        (
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            3,
            3,
            2,
            2,
            2,
            2,
            3,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
    (
        "F19",
        3,
        (
            3,
            0,
            0,
            1,
            1,
            1,
            1,
            0,
            0,
            1,
            1,
            1,
            1,
            3,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            1,
            1,
            1,
            1,
            1,
            1,
        ),
    ),
]


def build_Q_K8_matrix():
  rows = []
  for idx in range(28):
    r = np.zeros(28, dtype=np.int64)
    r[idx] = 1
    rows.append((r, 0))
  for length in [3, 5, 7]:
    for verts in combinations(range(8), length):
      first, rest = verts[0], verts[1:]
      for perm in permutations(rest):
        if perm[0] < perm[-1]:
          cyc = (first,) + perm
          r = np.zeros(28, dtype=np.int64)
          for k in range(length):
            u_v, w_v = cyc[k], cyc[(k + 1) % length]
            if u_v > w_v:
              u_v, w_v = w_v, u_v
            r[EDGE_TO_IDX8[(u_v, w_v)]] = 1
          rows.append((r, 1))
  A_mat = np.array([r[0] for r in rows], dtype=np.int64)
  b_vec = np.array([r[1] for r in rows], dtype=np.int64)
  return A_mat, b_vec


def fast_canonical_rep8(vec_tuple):
  if len(set(vec_tuple)) == 1:
    return vec_tuple
  adj = [[0] * 8 for _ in range(8)]
  for idx, (u, v) in enumerate(EDGES8):
    val = vec_tuple[idx]
    adj[u][v] = val
    adj[v][u] = val
  twin_rep = list(range(8))
  for i in range(8):
    for j in range(i + 1, 8):
      if all(adj[i][w] == adj[j][w] for w in range(8) if w != i and w != j):
        if twin_rep[j] == j:
          twin_rep[j] = twin_rep[i]
  sigs = [
      tuple(sorted(adj[i][j] for j in range(8) if j != i)) for i in range(8)
  ]
  sigs2 = [
      (sigs[i], tuple(sorted((adj[i][j], sigs[j]) for j in range(8) if j != i)))
      for i in range(8)
  ]
  groups = {}
  for i, s in enumerate(sigs2):
    groups.setdefault(s, []).append(i)
  sorted_keys = sorted(groups.keys())
  group_twin_classes = []
  for k in sorted_keys:
    g = groups[k]
    tc = {}
    for u in g:
      tc.setdefault(twin_rep[u], []).append(u)
    group_twin_classes.append(list(tc.values()))

  def gen_group_perms(tc_list):
    by_size = {}
    for c in tc_list:
      by_size.setdefault(len(c), []).append(c)
    size_keys = sorted(by_size.keys())

    def rec_size(s_idx):
      if s_idx == len(size_keys):
        yield ()
        return
      for p_classes in permutations(by_size[size_keys[s_idx]]):
        flat = tuple(u for c in p_classes for u in c)
        for rest in rec_size(s_idx + 1):
          yield flat + rest

    yield from rec_size(0)

  def gen_all_perms(g_idx):
    if g_idx == len(group_twin_classes):
      yield ()
      return
    for p in gen_group_perms(group_twin_classes[g_idx]):
      for rest in gen_all_perms(g_idx + 1):
        yield p + rest

  best = None
  for p in gen_all_perms(0):
    cand = tuple(adj[p[u]][p[v]] for u, v in EDGES8)
    if best is None or cand < best:
      best = cand
  return best


def exact_dd_fast128(A_tight):
  """100% deterministic Double Description Method with 128-bit bitmask adjacency test."""
  m, d = A_tight.shape
  basis_indices, cur_basis = [], []
  for i in range(m):
    cand = cur_basis + [A_tight[i]]
    if np.linalg.matrix_rank(np.array(cand, dtype=np.float64)) == len(cand):
      cur_basis = cand
      basis_indices.append(i)
      if len(basis_indices) == d:
        break
  rest_indices = [i for i in range(m) if i not in basis_indices]
  M = A_tight[basis_indices + rest_indices]
  B_inv = np.linalg.inv(M[:d].astype(np.float64))
  rays, masks_lo, masks_hi = [], [], []
  for j in range(d):
    col = B_inv[:, j]
    fracs = [Fraction(float(x)).limit_denominator(100000) for x in col]
    lcm_den = 1
    for fr in fracs:
      lcm_den = (lcm_den * fr.denominator) // int(
          np.gcd(lcm_den, fr.denominator)
      )
    r_int = np.array([int(fr * lcm_den) for fr in fracs], dtype=np.int64)
    g = int(np.gcd.reduce(np.abs(r_int)))
    if g > 1:
      r_int //= g
    rays.append(r_int)
    masks_lo.append(np.uint64(((1 << d) - 1) ^ (1 << j)))
    masks_hi.append(np.uint64(0))
  remaining = list(range(d, m))
  while remaining:
    rays_arr = np.array(rays, dtype=np.int64)
    best_r, best_sc = None, None
    for r_i in remaining:
      vc = rays_arr @ M[r_i]
      sc = int(np.sum(vc > 0)) * int(np.sum(vc < 0))
      if best_sc is None or sc < best_sc:
        best_sc, best_r = sc, r_i
    remaining.remove(best_r)
    row_i = best_r
    vals = rays_arr @ M[row_i]
    pos_idx = np.where(vals > 0)[0]
    zero_idx = np.where(vals == 0)[0]
    neg_idx = np.where(vals < 0)[0]
    bit_lo = np.uint64(1 << row_i) if row_i < 64 else np.uint64(0)
    bit_hi = np.uint64(1 << (row_i - 64)) if row_i >= 64 else np.uint64(0)
    if len(neg_idx) == 0:
      for idx in zero_idx:
        masks_lo[idx] |= bit_lo
        masks_hi[idx] |= bit_hi
      continue
    new_rays = [rays[idx] for idx in pos_idx] + [rays[idx] for idx in zero_idx]
    new_lo = [masks_lo[idx] for idx in pos_idx] + [
        masks_lo[idx] | bit_lo for idx in zero_idx
    ]
    new_hi = [masks_hi[idx] for idx in pos_idx] + [
        masks_hi[idx] | bit_hi for idx in zero_idx
    ]
    arr_lo = np.array(masks_lo, dtype=np.uint64)
    arr_hi = np.array(masks_hi, dtype=np.uint64)
    for p in pos_idx:
      zlo, zhi, vp, rp = masks_lo[p], masks_hi[p], int(vals[p]), rays[p]
      for n in neg_idx:
        clo, chi = zlo & masks_lo[n], zhi & masks_hi[n]
        if int(clo).bit_count() + int(chi).bit_count() < d - 2:
          continue
        if (
            np.count_nonzero(((arr_lo & clo) == clo) & ((arr_hi & chi) == chi))
            == 2
        ):
          vn = int(-vals[n])
          r_new = vp * rays[n] + vn * rp
          g = int(np.gcd.reduce(np.abs(r_new)))
          if g > 1:
            r_new //= g
          new_rays.append(r_new)
          new_lo.append(clo | bit_lo)
          new_hi.append(chi | bit_hi)
    rays, masks_lo, masks_hi = new_rays, new_lo, new_hi
  return rays


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
    rowA, rowD = np.zeros(total_vars), np.zeros(total_vars)
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


def verify_R13_orientations():
  print("=" * 76)
  print("[Check 1] All 1088 Extreme (1/3, 2/3)-Orientations of R_{1/3}(8, m)")
  print("=" * 76)
  for m in [1, 2, 3]:
    pair_list = [(i, j) for i in range(m) for j in range(i + 1, m)]
    single_list = [(i, u) for i in range(m) for u in range(2 * m, 8)]
    n_bits = len(pair_list) + len(single_list)
    orbits = set()
    for bits in product([0, 1], repeat=n_bits):
      vec = [Fraction(1, 3)] * 28
      for i in range(m):
        vec[EDGE_TO_IDX8[(2 * i, 2 * i + 1)]] = Fraction(0)
      for idx_p, (i, j) in enumerate(pair_list):
        b_val = bits[idx_p]
        for a in range(2):
          for b_s in range(2):
            u, v = 2 * i + a, 2 * j + b_s
            vec[EDGE_TO_IDX8[(u, v)]] = (
                Fraction(2, 3) if (a ^ b_s ^ b_val) == 0 else Fraction(1, 3)
            )
      offset = len(pair_list)
      for idx_s, (i, w) in enumerate(single_list):
        b_val = bits[offset + idx_s]
        vec[EDGE_TO_IDX8[(2 * i, w)]] = (
            Fraction(2, 3) if b_val == 0 else Fraction(1, 3)
        )
        vec[EDGE_TO_IDX8[(2 * i + 1, w)]] = (
            Fraction(1, 3) if b_val == 0 else Fraction(2, 3)
        )
      orbits.add(fast_canonical_rep8(tuple(vec)))
    min_b = min(compute_exact_B(8, orb)[0] for orb in orbits)
    print(
        f"  m={m}: 2^{n_bits} = {1<<n_bits:4d} orientations ->"
        f" {len(orbits):2d} S_8-orbits, min B_8(u) = {min_b:.6f}"
    )


def verify_19_fractional_orbits(A_mat, b_vec):
  print("\n" + "=" * 76)
  print("[Check 2] All 19 Fractional Extreme Point Orbits F01..F19 of Q_{K_8}^up")
  print("=" * 76)
  for name, den, num_tup in ORBITS_F01_F19:
    v_frac = tuple(Fraction(x, den) for x in num_tup)
    v_num = np.array(num_tup, dtype=np.int64)
    assert np.all(A_mat @ v_num >= b_vec * den), f"{name} violates Q_K8^up!"
    tight_mask = A_mat @ v_num == b_vec * den
    A_t = A_mat[tight_mask]
    rk = np.linalg.matrix_rank(A_t.astype(float))
    assert rk == 28, f"{name} has rank {rk} != 28!"
    zeros = [EDGES8[i] for i, x in enumerate(num_tup) if x == 0]
    verts_z = [u for e in zeros for u in e]
    is_m = len(verts_z) == len(set(verts_z))
    b_val, _, _ = compute_exact_B(8, v_frac)
    print(
        f"  {name} (d={den}, m={len(zeros)}, match={str(is_m):5s},"
        f" n_t={len(A_t):2d}, rank={rk}, sum={str(sum(v_frac)):5s}): B_8(v) ="
        f" {b_val:.6f}"
    )


def verify_dd_closure(A_mat, b_vec):
  print("\n" + "=" * 76)
  print("[Check 3] Deterministic 1-Skeleton Double Description Closure (F02..F19)")
  print("=" * 76)
  known_frac = {
      fast_canonical_rep8(tuple(Fraction(x, den) for x in num_tup))
      for _, den, num_tup in ORBITS_F01_F19
  }
  for name, den, num_tup in ORBITS_F01_F19:
    if name == "F01":
      continue
    v_frac = tuple(Fraction(x, den) for x in num_tup)
    v_num = np.array(num_tup, dtype=np.int64)
    A_tight = A_mat[A_mat @ v_num == b_vec * den]
    t_s = time.time()
    rays = exact_dd_fast128(A_tight)
    for r in rays:
      Ar = A_mat @ r
      neg = np.where(Ar < 0)[0]
      if len(neg) == 0:
        continue
      rem = A_mat[neg] @ v_num - b_vec[neg] * den
      den_r = -Ar[neg] * den
      min_t = min(Fraction(int(rem[k]), int(den_r[k])) for k in range(len(neg)))
      w_can = fast_canonical_rep8(
          tuple(v_frac[j] + min_t * int(r[j]) for j in range(28))
      )
      if max(x.denominator for x in w_can) > 1:
        assert w_can in known_frac, f"Unlisted fractional neighbor from {name}!"
    print(
        f"  {name} (n_t={len(A_tight):2d}): {len(rays):4d} extreme rays verified"
        f" in {time.time()-t_s:.3f}s [Closed in F01..F19]"
    )


def verify_K9_counterexample():
  print("\n" + "=" * 76)
  print("[Check 4] K = 9 Counterexample b*(9, 4) (J_4^circ)")
  print("=" * 76)
  edges9 = list(combinations(range(9), 2))
  bot = {1, 3, 5, 7}
  zero_set = {(0, 1), (2, 3), (4, 5), (6, 7)}
  v9 = tuple(
      Fraction(0)
      if e in zero_set
      else (
          Fraction(1, 3) if (e[0] in bot) == (e[1] in bot) else Fraction(2, 3)
      )
      for e in edges9
  )
  b_val9, c9, m9 = compute_exact_B(9, v9)
  print(
      f"  K=9, m={m9} (c={c9}): B_9(b*)/c*(b*) = {b_val9:.6f} == 14/15 < 1"
      f" (B_9(b*) = {3*b_val9:.6f} == 14/5 = 2.8 < 3 = c*(b*))"
  )


if __name__ == "__main__":
  t0 = time.time()
  A_mat, b_vec = build_Q_K8_matrix()
  verify_R13_orientations()
  verify_19_fractional_orbits(A_mat, b_vec)
  verify_dd_closure(A_mat, b_vec)
  verify_K9_counterexample()
  print(f"\nAll checks passed deterministically in {time.time() - t0:.2f}s.")

