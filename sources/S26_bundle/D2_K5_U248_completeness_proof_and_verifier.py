#!/usr/bin/env python3
"""D2, K=5: Exact 248-Term Closed Form (U248), U188 Counterexample Verifier,
12-Orbit Irredundancy Audit, and Complete 11-Polyhedron 1-Skeleton + Farkas Proof Suite.

Usage:
  python3 D2_K5_U248_completeness_proof_and_verifier.py             # Fast exact Fraction audit (~1s)
  python3 D2_K5_U248_completeness_proof_and_verifier.py --full-bfs  # Full 1-skeleton BFS + 335 Farkas proofs (~95s)
"""
from __future__ import annotations
import argparse
import itertools
import json
import time
from dataclasses import dataclass
from fractions import Fraction as F
import numpy as np

K = 5
EDGES = tuple(itertools.combinations(range(K), 2))
STATES = np.array([[(mask >> i) & 1 for i in range(K)] for mask in range(32)], dtype=int)
PAIR = np.array([STATES[:, i] * STATES[:, j] for i, j in EDGES]).T
MOMENT = np.column_stack([np.ones(32, dtype=int), STATES, PAIR])
I10 = np.hstack([np.zeros((10, 6), dtype=int), np.eye(10, dtype=int)])
A_BASE_INT = np.vstack([MOMENT, I10])  # 42 x 16 integer matrix


@dataclass(frozen=True)
class Certificate:
    label: str
    a: F
    beta: tuple[F, ...]
    gamma: tuple[F, ...]
    weights: tuple[int, ...]

    def coefficients(self) -> tuple[F, ...]:
        return (self.a,) + self.beta + self.gamma

    def bound(self, q, b) -> F:
        return (
            self.a
            + sum(v * F(x) for v, x in zip(self.beta, q))
            + sum(v * F(x) for v, x in zip(self.gamma, b))
        )

    def check(self) -> None:
        assert len(self.beta) == 5 and len(self.gamma) == 10
        assert all(z >= 0 for z in self.gamma)
        assert sum(self.weights) % 2 == 1 and sum(self.weights) > 0
        for m in range(32):
            x = STATES[m]
            lhs = (
                self.a
                + sum(c * int(t) for c, t in zip(self.beta, x))
                + sum(c * int(PAIR[m, e]) for e, c in enumerate(self.gamma))
            )
            rhs = int(2 * sum(int(v) * int(t) for v, t in zip(self.weights, x)) > sum(self.weights))
            assert lhs >= rhs, (self.label, m, lhs, rhs)


def apply_perm(v_rat: tuple[F, ...], p: tuple[int, ...]) -> tuple[F, ...]:
    beta_new = tuple(v_rat[1 + p[i]] for i in range(5))
    gamma_new = tuple(v_rat[6 + EDGES.index(tuple(sorted((p[i], p[j]))))] for i, j in EDGES)
    return (v_rat[0],) + beta_new + gamma_new


def odd_cycles():
    cc = []
    for vertices in itertools.combinations(range(K), 3):
        cc.append(tuple(EDGES.index(e) for e in itertools.combinations(vertices, 2)))
    for rest in itertools.permutations(range(1, K)):
        cyc = (0,) + rest
        if cyc[1] > cyc[-1]:
            continue
        cc.append(tuple(EDGES.index(tuple(sorted((cyc[i], cyc[(i + 1) % K])))) for i in range(K)))
    assert len(cc) == 22
    return cc


def weighted_bp(weights, label):
    w = tuple(weights)
    total = sum(w)
    r = (total + 1) // 2
    a = F((r - 2) * (r - 1), 2)
    beta = tuple(F(ww * (ww - 2 * r + 3), 2) for ww in w)
    gamma = tuple(F(w[i] * w[j]) for i, j in EDGES)
    return Certificate(label, a, beta, gamma, w)


def rooted(i):
    w = tuple(3 if j == i else 1 for j in range(K))
    a = F(1, 3)
    beta = (F(-1, 3),) * K
    gamma = tuple(F(4 if i in (j, k) else 1, 3) for j, k in EDGES)
    return Certificate(f"rooted:{i}", a, beta, gamma, w)


def corrected_core(weights):
    w = tuple(weights)
    root = weighted_bp(w, "unused")
    H = {i for i, wi in enumerate(w) if wi >= 2}
    beta = tuple(root.beta[i] + (1 if i in H else 0) for i in range(K))
    gamma = tuple(root.gamma[e] - (1 if i in H and j in H else 0) for e, (i, j) in enumerate(EDGES))
    return Certificate("corrected:" + "".join(map(str, w)), root.a - 1, beta, gamma, w)


def make_u98() -> list[Certificate]:
    out = []
    for i in range(K):
        w = tuple(int(j == i) for j in range(K))
        out.append(Certificate(f"dict:{i}", F(0), tuple(F(v) for v in w), (F(0),) * 10, w))
    for ix, C_edges in enumerate(odd_cycles()):
        vertices = {i for e in C_edges for i in EDGES[e]}
        w = tuple(int(i in vertices) for i in range(K))
        gamma = tuple(F(int(e in C_edges)) for e in range(10))
        out.append(Certificate(f"cycle:{ix}", F(0), (F(0),) * 5, gamma, w))
    for canonical in ((3, 1, 1, 1, 1), (2, 1, 1, 1, 0), (2, 2, 1, 1, 1), (1, 1, 1, 1, 1)):
        for w in sorted(set(itertools.permutations(canonical))):
            out.append(weighted_bp(w, "bp:" + "".join(map(str, w))))
    for i in range(K):
        out.append(rooted(i))
    for w in sorted(set(itertools.permutations((3, 2, 2, 1, 1)))):
        out.append(corrected_core(w))
    assert len(out) == 98
    return out


def make_anchored30() -> list[Certificate]:
    out = []
    for r in range(5):
        others = [i for i in range(5) if i != r]
        for s, t in itertools.combinations(others, 2):
            weights = tuple(3 if i == r else 1 if i in (s, t) else 2 for i in range(5))
            beta = tuple(F(0) if i == r else F(-1) for i in range(5))
            gamma = tuple(
                F(2) if set((i, j)) == set((s, t))
                else F(0) if r in (i, j) and (s in (i, j) or t in (i, j))
                else F(1)
                for i, j in EDGES
            )
            out.append(Certificate(f"anchored2:{r}:{s}-{t}", F(1), beta, gamma, weights))
    assert len(out) == 30
    return out


def make_twoport60() -> list[Certificate]:
    result = []
    for h in range(5):
        for medium in itertools.combinations([i for i in range(5) if i != h], 2):
            light = [i for i in range(5) if i != h and i not in medium]
            H = (h, *light)
            for assignment in itertools.permutations(light):
                cross = {
                    tuple(sorted((medium[0], h))),
                    tuple(sorted((medium[1], h))),
                    tuple(sorted((medium[0], assignment[0]))),
                    tuple(sorted((medium[1], assignment[1]))),
                }
                gamma = tuple(F(int(set(ij) <= set(H) or ij in cross)) for ij in EDGES)
                beta = tuple(F(-1) if i in H else F(0) for i in range(5))
                weights = tuple(3 if i == h else 2 if i in medium else 1 for i in range(5))
                result.append(
                    Certificate(f"twoport:root={h}:medium={medium}:map={assignment}", F(1), beta, gamma, weights)
                )
    assert len(result) == 60
    return result


def make_coupled_pair60() -> list[Certificate]:
    """The newly discovered 60-certificate S5-orbit of 9-edge quadratic majorants."""
    v_seed = (
        F(2),
        F(-2), F(-1), F(-2), F(-1), F(-2),
        F(2), F(3), F(1), F(2), F(1), F(0), F(1), F(2), F(2), F(1),
    )
    orb = sorted({apply_perm(v_seed, p) for p in itertools.permutations(range(5))})
    assert len(orb) == 60
    out = []
    for idx, coeff in enumerate(orb):
        beta = coeff[1:6]
        gamma = coeff[6:]
        e3 = next(e for e, g in enumerate(gamma) if g == 3)
        a1, a2 = EDGES[e3]
        # l1 is the node with beta == -1 and gamma_{a1, l1} == 2; majorizes both (2,2,1,1,1) on A={a1,a2} and (3,2,2,1,1) on {a1:3, a2:2, l1:2}
        l1 = next(i for i in range(5) if beta[i] == -1 and gamma[EDGES.index(tuple(sorted((a1, i))))] == 2)
        weights = tuple(3 if i == a1 else 2 if i in (a2, l1) else 1 for i in range(5))
        cert = Certificate(f"coupled_pair:{idx}", coeff[0], beta, gamma, weights)
        out.append(cert)
    return out


def make_u188() -> list[Certificate]:
    return make_u98() + make_anchored30() + make_twoport60()


def make_u248() -> list[Certificate]:
    certs = make_u188() + make_coupled_pair60()
    assert len(certs) == 248
    assert len({c.coefficients() for c in certs}) == 248
    for c in certs:
        c.check()
    return certs


def bound_u248(q, b) -> tuple[F, str]:
    return min((c.bound(q, b), c.label) for c in make_u248())


def verify_u188_counterexamples() -> dict:
    """Verify exact rational SAT counterexamples disproving U188 on both 8-edge branches and on 10-edge exact moments."""
    u188 = make_u188()
    v_rat = (
        F(2),
        F(-2), F(-1), F(-2), F(-1), F(-2),
        F(2), F(3), F(1), F(2), F(1), F(0), F(1), F(2), F(2), F(1),
    )
    new_to_old = [1, 3, 0, 2, 4]
    pi0 = [F(1, 6) if m in (1, 6, 9, 17, 20, 26) else F(0) for m in range(32)]
    p0 = [sum(pi0[m] * int(PAIR[m, e]) for m in range(32)) for e in range(10)]
    b0 = list(p0)
    b0[5] = F(1, 3)

    delta = F(1, 100)
    anchor_pi = [F(1, 5) ** int(sum(STATES[m])) * F(4, 5) ** int(5 - sum(STATES[m])) for m in range(32)]
    pi_int = [(1 - delta) * pi0[m] + delta * anchor_pi[m] for m in range(32)]
    p_int = [sum(pi_int[m] * int(PAIR[m, e]) for m in range(32)) for e in range(10)]
    b_int_old = [p_int[e] if e != 5 else (1 - delta) * b0[5] + delta * F(1, 25) for e in range(10)]

    pi_9sat = [F(0)] * 32
    for m_new in range(32):
        m_old = sum((1 << new_to_old[i]) for i in range(5) if (m_new >> i) & 1)
        pi_9sat[m_new] = pi_int[m_old]
    b_9sat = [b_int_old[EDGES.index(tuple(sorted((new_to_old[i], new_to_old[j]))))] for i, j in EDGES]
    gamma_9sat = [v_rat[6 + EDGES.index(tuple(sorted((new_to_old[i], new_to_old[j]))))] for i, j in EDGES]
    beta_9sat = [v_rat[1 + new_to_old[i]] for i in range(5)]
    Q_9sat = (v_rat[0], *beta_9sat, *gamma_9sat)
    q_9sat = [sum(pi_9sat[m] * int(STATES[m, i]) for m in range(32)) for i in range(5)]
    p_9sat = [sum(pi_9sat[m] * int(PAIR[m, e]) for m in range(32)) for e in range(10)]

    assert all(p > 0 for p in pi_9sat) and sum(pi_9sat) == 1
    assert all(F(0) < q < F(1, 2) for q in q_9sat) and sum(q_9sat) - 3 * min(q_9sat) < 1
    assert p_9sat[0] < b_9sat[0] < min(q_9sat[0], q_9sat[1])
    assert all(p_9sat[e] == b_9sat[e] for e in range(1, 10))
    assert all(Q_9sat[6 + e] > 0 for e in range(1, 10)) and Q_9sat[6] == 0

    u188_9sat = min(c.bound(q_9sat, b_9sat) for c in u188)
    lv_9sat = Q_9sat[0] + sum(Q_9sat[1 + i] * q_9sat[i] for i in range(5)) + sum(Q_9sat[6 + e] * b_9sat[e] for e in range(10))
    gap_9sat = u188_9sat - lv_9sat
    assert u188_9sat == F(831, 5000) and lv_9sat == F(1, 100) and gap_9sat == F(781, 5000)

    c_corr = corrected_core((3, 2, 2, 1, 1)).coefficients()
    eps = F(1, 1000)
    Q_10sat = tuple((1 - eps) * v_rat[k] + eps * c_corr[k] for k in range(16))
    q_10sat = [sum(pi_int[m] * int(STATES[m, i]) for m in range(32)) for i in range(5)]
    b_10sat = list(p_int)
    u188_10sat = min(c.bound(q_10sat, b_10sat) for c in u188)
    lv_10sat = Q_10sat[0] + sum(Q_10sat[1 + i] * q_10sat[i] for i in range(5)) + sum(Q_10sat[6 + e] * b_10sat[e] for e in range(10))
    gap_10sat = u188_10sat - lv_10sat
    assert u188_10sat == F(17, 1250) and lv_10sat == F(25033, 2500000) and gap_10sat == F(8967, 2500000)

    return {
        "9sat_both_branches": {"U188": str(u188_9sat), "L_Q": str(lv_9sat), "gap": str(gap_9sat)},
        "10sat_exact_moments": {"U188": str(u188_10sat), "L_Q": str(lv_10sat), "gap": str(gap_10sat)},
    }


SPARSE_12_ORBIT_WITNESSES = [
    (0, 'dict:0', 5, {0: '1/2', 6: '1/4', 24: '1/4'}, {4: '1/4', 5: '1/4', 6: '1/4', 7: '1/4', 8: '1/4', 9: '1/4'}),
    (5, 'cycle:0', 10, {1: '1/6', 2: '1/6', 9: '1/6', 12: '1/6', 18: '1/6', 20: '1/6'}, {2: '1/6', 3: '1/6', 5: '1/6', 6: '1/6', 7: '1/6', 8: '1/6'}),
    (15, 'cycle:10', 12, {5: '1/5', 9: '1/5', 10: '1/5', 18: '1/5', 20: '1/5'}, {1: '2/5', 2: '2/5', 5: '2/5', 6: '2/5', 8: '2/5'}),
    (27, 'bp:11113', 5, {7: '3/20', 11: '3/20', 13: '3/20', 14: '3/20', 16: '2/5'}, {0: '3/10', 1: '3/10', 2: '3/10', 4: '3/10', 5: '3/10', 7: '3/10'}),
    (32, 'bp:01112', 20, {6: '2/9', 10: '2/9', 12: '1/9', 13: '1/9', 16: '1/9', 17: '2/9'}, {0: '1/9', 1: '1/9', 2: '1/9', 3: '2/9', 4: '2/9', 5: '2/9', 7: '2/9'}),
    (52, 'bp:11122', 10, {7: '1/5', 8: '1/10', 9: '1/10', 10: '1/10', 12: '1/10', 16: '1/10', 17: '1/10', 18: '1/10', 20: '1/10'}, {0: '1/5', 1: '1/5', 2: '1/10', 3: '1/10', 4: '1/5', 5: '1/10', 6: '1/10', 7: '1/10', 8: '1/10'}),
    (62, 'bp:11111', 1, {1: '1/20', 2: '1/20', 3: '3/40', 4: '1/20', 5: '3/40', 6: '3/40', 8: '1/20', 9: '3/40', 10: '3/40', 12: '3/40', 16: '1/20', 17: '3/40', 18: '3/40', 20: '3/40', 24: '3/40'}, {0: '3/40', 1: '3/40', 2: '3/40', 3: '3/40', 4: '3/40', 5: '3/40', 6: '3/40', 7: '3/40', 8: '3/40', 9: '3/40'}),
    (63, 'rooted:0', 5, {1: '5/26', 2: '3/26', 6: '3/26', 8: '3/26', 12: '3/26', 18: '3/26', 20: '3/26', 24: '3/26'}, {4: '3/26', 5: '3/26', 6: '3/26', 7: '3/26', 8: '3/26', 9: '3/26'}),
    (68, 'corrected:11223', 30, {7: '1/6', 11: '1/6', 12: '1/4', 16: '1/12', 17: '1/6', 18: '1/6'}, {0: '1/3', 1: '1/6', 2: '1/6', 3: '1/6', 4: '1/6', 5: '1/6', 6: '1/6', 7: '1/4'}),
    (98, 'anchored2:0:1-2', 30, {3: '2/9', 5: '1/9', 10: '1/9', 12: '1/9', 18: '1/9', 20: '1/9', 24: '2/9'}, {0: '1/3', 1: '2/9', 5: '1/9', 6: '1/9', 7: '1/9', 8: '1/9', 9: '2/9'}),
    (128, 'twoport:root=0:medium=(1, 2):map=(3, 4)', 60, {1: '1/14', 9: '1/7', 12: '2/7', 17: '1/7', 18: '2/7', 30: '1/14'}, {2: '1/7', 3: '1/7', 4: '3/14', 5: '1/14', 6: '5/14', 7: '5/14', 8: '1/14', 9: '1/14'}),
    (188, 'coupled_pair:0', 60, {2: '1/6', 3: '1/6', 5: '1/6', 10: '1/6', 20: '1/6', 25: '1/6'}, {0: '1/6', 1: '1/6', 2: '1/6', 3: '1/6', 5: '1/6', 8: '1/6', 9: '1/3'}),
]


def verify_12_orbit_irredundancy() -> list[dict]:
    """Verify in pure Fraction arithmetic that all 12 S5-orbits (all 248 terms of U248) are strictly exposed on int(R)."""
    u248 = make_u248()
    delta = F(1, 1000)
    anchor_pi = [F(1, 5) ** int(sum(STATES[m])) * F(4, 5) ** int(5 - sum(STATES[m])) for m in range(32)]
    anchor_b = [F(1, 25) + F(1, 100) for _ in range(10)]
    results = []
    total_sz = 0
    for idx, label, sz, pi_nz, b_nz in SPARSE_12_ORBIT_WITNESSES:
        assert u248[idx].label == label
        vv_pi = [F(pi_nz.get(m, '0')) for m in range(32)]
        vv_b = [F(b_nz.get(e, '0')) for e in range(10)]
        pm = [(1 - delta) * vv_pi[m] + delta * anchor_pi[m] for m in range(32)]
        qm = [sum(pm[m] * int(STATES[m, k]) for m in range(32)) for k in range(5)]
        bm = [(1 - delta) * vv_b[e] + delta * anchor_b[e] for e in range(10)]
        assert sum(pm) == 1 and all(p > 0 for p in pm)
        assert all(F(0) < q < F(1, 2) for q in qm) and all(sum(qm) - 3 * q < 1 for q in qm)
        assert all(sum(pm[m] * int(PAIR[m, e]) for m in range(32)) < bm[e] < min(qm[u], qm[v]) for e, (u, v) in enumerate(EDGES))
        val_this = u248[idx].bound(qm, bm)
        min_other = min(u248[j].bound(qm, bm) for j in range(248) if j != idx)
        gap = min_other - val_this
        assert gap > 0, (idx, label, gap)
        total_sz += sz
        results.append({"index": idx, "label": label, "orbit_size": sz, "strict_gap": str(gap)})
    assert total_sz == 248
    return results


def tangent_extreme_rays_gen(A_mat_int: np.ndarray, act_indices: tuple[int, ...]) -> list[np.ndarray]:
    M = A_mat_int[list(act_indices)]
    m, d = M.shape
    basis_rows, rem_rows = [], []
    cur_mat = np.empty((0, d), dtype=float)
    for idx in range(m):
        cand = np.vstack([cur_mat, M[idx]])
        if np.linalg.matrix_rank(cand) > len(basis_rows):
            basis_rows.append(idx)
            cur_mat = cand
        else:
            rem_rows.append(idx)
    B_inv = np.linalg.inv(M[basis_rows].astype(float))
    rays = []
    for col in range(16):
        r = B_inv[:, col]
        r /= np.max(np.abs(r))
        zmask = sum((1 << basis_rows[i]) for i in range(16) if i != col)
        rays.append((r, zmask))

    rem = list(rem_rows)
    while rem:
        best_r, best_score = None, None
        for r_idx in rem:
            row = M[r_idx].astype(float)
            p_cnt = sum(1 for r, _ in rays if row @ r > 1e-9)
            n_cnt = sum(1 for r, _ in rays if row @ r < -1e-9)
            score = p_cnt * n_cnt
            if best_score is None or score < best_score:
                best_score = score
                best_r = r_idx
        rem.remove(best_r)
        row = M[best_r].astype(float)
        pos, neg, zer = [], [], []
        bit = 1 << best_r
        for r, zmask in rays:
            val = float(row @ r)
            if abs(val) < 1e-9:
                zer.append((r, zmask | bit))
            elif val > 0:
                pos.append((r, zmask, val))
            else:
                neg.append((r, zmask, -val))
        new_rays = [(r, zmask) for r, zmask, _ in pos] + list(zer)
        all_zmasks = [z for _, z in rays]
        for rp, zp, vp in pos:
            for rn, zn, vn in neg:
                common = zp & zn
                if common.bit_count() < 14:
                    continue
                if sum(1 for zw in all_zmasks if (zw & common) == common) == 2:
                    r_new = vn * rp + vp * rn
                    r_new /= np.max(np.abs(r_new))
                    new_rays.append((r_new, common | bit))
        rays = new_rays
    return [r for r, _ in rays]


def run_full_11polyhedra_bfs_and_farkas_audit() -> dict:
    """Run untruncated 1-seed 1-skeleton BFS on all 11 canonical dual polyhedra (7 monotone + 4 drop-0)
    and verify exact rational Farkas dominance certificates for all 335 non-U248 S5-orbits."""
    from scipy.optimize import linprog

    u248 = [c.coefficients() for c in make_u248()]
    u248_set = set(u248)
    coeff248 = np.array([[float(t) for t in c] for c in u248])
    Bmom = MOMENT[:, 1:6].astype(float)
    coefmat248 = np.hstack((coeff248[:, 1:6] @ Bmom.T, coeff248[:, 6:], np.zeros((248, 1))))

    base, baserhs = [], []
    def add_row(pi_r, b_r, r_val):
        base.append(tuple([F(t) for t in pi_r + b_r] + [F(0)]))
        baserhs.append(F(r_val))
    zero32, zero10 = [0] * 32, [0] * 10
    for e, (i, j) in enumerate(EDGES):
        z = zero10.copy(); z[e] = -1; add_row([int(PAIR[m, e]) for m in range(32)], z, 0)
        for k in (i, j):
            z = zero10.copy(); z[e] = 1; add_row([-int(STATES[m, k]) for m in range(32)], z, 0)
    for k in range(5):
        add_row([int(STATES[m, k]) for m in range(32)], zero10, F(1, 2))
        add_row([int(sum(STATES[m]) - 3 * STATES[m, k]) for m in range(32)], zero10, 1)
    for e in range(10):
        z = zero10.copy(); z[e] = -1; add_row(zero32, z, 0)
    eq = [[F(1)] * 32 + [F(0)] * 11]
    eqrhs = [F(1)]
    basef = np.array(base, dtype=float)
    b0 = np.array(baserhs, dtype=float)
    eqf = np.array(eq, dtype=float)
    eq0 = np.array(eqrhs, dtype=float)

    def rat_vec(x, maxden=10000):
        return tuple(F(float(v)).limit_denominator(maxden) for v in x)

    def bfs_polyhedron(A_sys_int, b_sys_int, auts, seeds):
        A_sys = A_sys_int.astype(float)
        nrows = len(b_sys_int)
        def get_act(v_r):
            act = []
            for idx in range(nrows):
                lhs = sum(int(A_sys_int[idx, k]) * v_r[k] for k in range(16))
                rhs = F(int(b_sys_int[idx]))
                assert lhs >= rhs
                if lhs == rhs:
                    act.append(idx)
            return tuple(act)
        valid_seeds = []
        for s in seeds:
            v_fl = np.array([float(x) for x in s])
            if np.min(A_sys @ v_fl - b_sys_int) >= -1e-9:
                act = get_act(s)
                if np.linalg.matrix_rank(A_sys[list(act)]) == 16:
                    valid_seeds.append((len(act), min(apply_perm(s, p) for p in auts)))
        valid_seeds.sort()
        seed_rep = valid_seeds[0][1]
        discovered_reps = {seed_rep}
        bfs_q = [seed_rep]
        head = 0
        while head < len(bfs_q):
            v_r = bfs_q[head]
            head += 1
            act = get_act(v_r)
            if len(act) == nrows:
                continue
            rays = tangent_extreme_rays_gen(A_sys_int, act)
            slacks_rat = [sum(int(A_sys_int[idx, k]) * v_r[k] for k in range(16)) - F(int(b_sys_int[idx])) for idx in range(nrows)]
            for r in rays:
                r_rat = rat_vec(r, 10000)
                Ad = [sum(int(A_sys_int[idx, k]) * r_rat[k] for k in range(16)) for idx in range(nrows)]
                if any(Ad[idx] < 0 for idx in act):
                    continue
                t_star = None
                for idx in range(nrows):
                    if Ad[idx] < 0 and slacks_rat[idx] > 0:
                        cand_t = slacks_rat[idx] / (-Ad[idx])
                        if t_star is None or cand_t < t_star:
                            t_star = cand_t
                if t_star is not None and t_star > 0:
                    v_next = tuple(vk + t_star * rk for vk, rk in zip(v_r, r_rat))
                    rep = min(apply_perm(v_next, p) for p in auts)
                    if rep not in discovered_reps:
                        act_next = get_act(rep)
                        assert np.linalg.matrix_rank(A_sys[list(act_next)]) == 16
                        discovered_reps.add(rep)
                        bfs_q.append(rep)
        return discovered_reps

    CANONICAL_WEIGHTS = [
        (1, 0, 0, 0, 0), (1, 1, 1, 0, 0), (1, 1, 1, 1, 1),
        (2, 1, 1, 1, 0), (2, 2, 1, 1, 1), (3, 1, 1, 1, 1), (3, 2, 2, 1, 1),
    ]
    perms = list(itertools.permutations(range(5)))
    seen_s5 = set()
    all_non_u248 = []

    for cw in CANONICAL_WEIGHTS:
        auts = [p for p in perms if tuple(cw[p[i]] for i in range(5)) == cw]
        ft = np.array([int(2 * sum(cw[i] * int(x[i]) for i in range(5)) > sum(cw)) for x in STATES], dtype=int)
        b_all = np.r_[ft, np.zeros(10, dtype=int)]
        reps = bfs_polyhedron(A_BASE_INT, b_all, auts, u248)
        for r in sorted(reps):
            if r not in u248_set:
                sc = min(apply_perm(r, p) for p in perms)
                if sc not in seen_s5:
                    seen_s5.add(sc)
                    all_non_u248.append((f"mon:{cw}", r))

    FU_CASES = [
        ("const_1", (0, 0, 0, 0), lambda x: 1),
        ("dict_4",  (1, 0, 0, 0), lambda x: int(x[1] == 1)),
        ("maj3_4",  (1, 1, 1, 0), lambda x: int(x[1] + x[2] + x[3] >= 2)),
        ("apex4_4", (2, 1, 1, 1), lambda x: int(2 * x[1] + x[2] + x[3] + x[4] >= 3)),
    ]
    rng = np.random.default_rng(42)
    for name, w4, fn in FU_CASES:
        auts = [p for p in perms if p[0] == 0 and tuple(w4[p[i] - 1] for i in range(1, 5)) == w4]
        ft = np.array([fn(x) if x[0] == 0 else 0 for x in STATES], dtype=int)
        pairB_rows = [MOMENT[m | 1] + MOMENT[(m | 1) ^ 31] for m in range(32) if (m & 1) == 0 and fn(STATES[m]) == 1]
        if name == "const_1":
            A_sys_int = A_BASE_INT
            b_sys_int = np.r_[ft, np.zeros(10, dtype=int)]
        else:
            A_sys_int = np.vstack([A_BASE_INT, np.array(pairB_rows, dtype=int)])
            b_sys_int = np.r_[ft, np.zeros(10, dtype=int), np.ones(len(pairB_rows), dtype=int)]
        seeds = []
        for _ in range(25):
            q = rng.uniform(0.05, 0.45, size=5)
            b = np.array([q[i] * q[j] for i, j in EDGES])
            res = linprog(np.r_[1.0, q, b], A_ub=-A_sys_int.astype(float), b_ub=-b_sys_int.astype(float), bounds=[(None, None)] * 16, method="highs")
            if res.success:
                seeds.append(rat_vec(res.x, 10000))
        reps = bfs_polyhedron(A_sys_int, b_sys_int, auts, seeds)
        for r in sorted(reps):
            if r not in u248_set:
                sc = min(apply_perm(r, p) for p in perms)
                if sc not in seen_s5:
                    seen_s5.add(sc)
                    all_non_u248.append((f"drop0:{name}", r))

    max_denom = 1
    for tag, v_r in all_non_u248:
        v_fl = np.array([float(x) for x in v_r])
        v_row = np.r_[v_fl[1:6] @ Bmom.T, v_fl[6:], 0.0]
        rj = v_row - coefmat248
        rj[:, -1] = 1.0
        bj = coeff248[:, 0] - v_fl[0]
        sol = linprog(
            np.r_[np.zeros(42), -1.0],
            A_ub=np.vstack((basef, rj)),
            b_ub=np.r_[b0, bj],
            A_eq=eqf,
            b_eq=eq0,
            bounds=[(0, None)] * 42 + [(None, None)],
            method="highs",
        )
        assert sol.success and -sol.fun <= 1e-8
        mu = [(ix, F(float(-v)).limit_denominator(100000)) for ix, v in enumerate(sol.ineqlin.marginals) if v < -1e-9]
        nu = [F(float(-v)).limit_denominator(100000) for v in sol.eqlin.marginals]
        accum = [F(0)] * 43
        cost = F(0)
        for ix, m in mu:
            assert m >= 0
            if ix < 50:
                row, r_val = base[ix], baserhs[ix]
            else:
                cj = u248[ix - 50]
                delta = [v_r[1 + k] - cj[1 + k] for k in range(5)]
                bb = [v_r[6 + e] - cj[6 + e] for e in range(10)]
                r_val = cj[0] - v_r[0]
                row = tuple([sum(delta[k] * int(STATES[m_idx, k]) for k in range(5)) for m_idx in range(32)] + bb + [F(1)])
            cost += m * r_val
            for k, val in enumerate(row):
                accum[k] += m * val
            max_denom = max(max_denom, m.denominator)
        for j, v_nu in enumerate(nu):
            cost += v_nu * eqrhs[j]
            for k, val in enumerate(eq[j]):
                accum[k] += v_nu * val
            max_denom = max(max_denom, v_nu.denominator)
        assert accum[-1] == 1 and min(accum[:-1]) >= 0 and cost <= 0

    return {"certified_non_u248_s5_orbits": len(all_non_u248), "max_farkas_denominator": max_denom}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full-bfs", action="store_true", help="Run full 11-polyhedron 1-skeleton BFS and 335 Farkas dominance audits")
    args = ap.parse_args()

    t0 = time.time()
    u248 = make_u248()
    print(f"[PASS] Pointwise check of all {len(u248)} certificates in U248 across all 32 states in {time.time() - t0:.2f}s.")

    sat_res = verify_u188_counterexamples()
    print(f"[PASS] Verified exact rational SAT counterexamples disproving U188: {json.dumps(sat_res)}")

    irred_res = verify_12_orbit_irredundancy()
    print(f"[PASS] Verified strict-interior exposure (irredundancy) of all {len(irred_res)} S5-orbits (all 248 terms of U248).")

    if args.full_bfs:
        t1 = time.time()
        bfs_res = run_full_11polyhedra_bfs_and_farkas_audit()
        print(f"[PASS] Full 11-polyhedron 1-skeleton BFS + exact rational Farkas audit completed in {time.time() - t1:.2f}s: {bfs_res}")


if __name__ == "__main__":
    main()
