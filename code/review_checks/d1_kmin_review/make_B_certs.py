"""Generate exact rational certificates of B_graph(u) >= 1 on K_8 (format checked by verify_B_certs.py).

Cutting-plane branch and bound.  Each node carries the path constraints y(T_j) >= 1.  The node LP
    min u.y  s.t.  y(T_j) >= 1, y >= 0
is solved EXACTLY in Fractions through its dual  max sum(lam)  s.t.  sum_j lam_j 1[T_j] <= u, lam >= 0
(dense tableau simplex, Bland's rule).  If the exact value is >= 1 the optimal dual lam is the leaf;
otherwise the exact primal optimum y* (read from the slack reduced costs) violates some cut S, i.e.
max(y*(E(S)), y*(E(V\\S))) < 1, and we branch on it.  Among violated cuts we pick by strong branching
(HiGHS float LP values of the two children; maximise the weaker child), tie-breaking towards the most
violated cut and towards cuts that separate the endpoints of every zero edge.

Usage: python3 -I make_B_certs.py <S27 orbit file> <out_dir> [orbit names...]   (skips existing files)
"""
import sys, os, ast, json, itertools, time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog

E = list(itertools.combinations(range(8), 2))
NE = len(E)
MAXCAND = int(os.environ.get('MAXCAND', '1000'))


def load_orbits(path):
    src = open(path).read()
    a = src.index('ORBITS_F01_F19 = [')
    b = src.index('\n]\n', a) + 2
    return {name: [F(x, d) for x in t] for name, d, t in ast.literal_eval(src[a + len('ORBITS_F01_F19 = '):b])}


def side(S, which):
    inS = [(S >> v) & 1 for v in range(8)]
    want = 1 if which == 'A' else 0
    return tuple(int(inS[i] == want and inS[j] == want) for i, j in E)


CUTS = [S for S in range(1, 255) if not (S >> 7) & 1]  # S and its complement give the same branching
SIDES = {S: (side(S, 'A'), side(S, 'D')) for S in CUTS}
SIDES_NP = {S: (np.array(a, float), np.array(d, float)) for S, (a, d) in SIDES.items()}


def exact_lp(T, u):
    """max sum lam s.t. sum_j lam_j T_j[e] <= u_e, lam >= 0.  Returns (value, lam, y)."""
    k = len(T)
    m = NE
    n = k + m
    # rows: constraint e: sum_j T_j[e] lam_j + s_e = u_e
    tab = [[F(T[j][e]) for j in range(k)] + [F(int(i == e)) for i in range(m)] + [u[e]] for e in range(m)]
    obj = [F(-1)] * k + [F(0)] * m + [F(0)]  # reduced costs (z - c), maximise
    basis = [k + e for e in range(m)]
    while True:
        enter = next((c for c in range(n) if obj[c] < 0), None)  # Bland
        if enter is None:
            break
        best = None
        for r in range(m):
            a = tab[r][enter]
            if a > 0:
                ratio = tab[r][n] / a
                if best is None or ratio < best[0] or (ratio == best[0] and basis[r] < basis[best[1]]):
                    best = (ratio, r)
        if best is None:
            raise RuntimeError('unbounded dual (cannot happen: every T_j is a nonempty edge set... unless u=0 on it)')
        r = best[1]
        piv = tab[r][enter]
        row = [x / piv for x in tab[r]]
        tab[r] = row
        for rr in range(m):
            if rr != r and tab[rr][enter] != 0:
                f = tab[rr][enter]
                tab[rr] = [x - f * y for x, y in zip(tab[rr], row)]
        if obj[enter] != 0:
            f = obj[enter]
            obj = [x - f * y for x, y in zip(obj, row)]
        basis[r] = enter
    lam = [F(0)] * k
    for r, bv in enumerate(basis):
        if bv < k:
            lam[bv] = tab[r][n]
    y = [obj[k + e] for e in range(m)]
    return obj[n], lam, y


def float_lp(T, u_np):
    if not T:
        return 0.0
    A = -np.array(T, float)
    r = linprog(u_np, A_ub=A, b_ub=-np.ones(len(T)), bounds=(0, None), method='highs')
    if r.status != 0:
        return float('inf')
    return r.fun


def choose_cut(T, y, u, u_np, zero_edges):
    pathset = set(T)
    cand = []
    for S in CUTS:
        a, d = SIDES[S]
        if a in pathset or d in pathset:
            continue
        ya = sum(y[e] for e in range(NE) if a[e])
        yd = sum(y[e] for e in range(NE) if d[e])
        mx = max(ya, yd)
        if mx < 1:
            sep = all(not a[e] and not d[e] for e in zero_edges)
            cand.append((S, 1 - mx, sep))
    assert cand, 'no violated cut although LP value < 1 (would contradict B_graph(u) >= 1 ... check)'
    cand.sort(key=lambda c: (c[2], c[1]), reverse=True)
    best = None
    for S, viol, sep in cand[:MAXCAND]:
        a, d = SIDES[S]
        va = float_lp(T + [a], u_np)
        vd = float_lp(T + [d], u_np)
        key = (round(min(va, vd), 9), round(va + vd, 9), sep, float(viol))
        if best is None or key > best[0]:
            best = (key, S)
    return best[1]


def build(T, u, u_np, zero_edges, stats, depth=0):
    stats['nodes'] += 1
    if T and not any(T[-1]):  # empty side (singleton cut): y(empty) >= 1 is impossible, lam = 1 on it
        stats['leaves'] += 1
        stats['depth'] = max(stats['depth'], len(T))
        return {'leaf': [[len(T) - 1, '1']]}
    val, lam, y = exact_lp(T, u)
    if val >= 1:
        stats['leaves'] += 1
        stats['depth'] = max(stats['depth'], len(T))
        return {'leaf': [[j, str(l)] for j, l in enumerate(lam) if l != 0]}
    if depth > 60:
        raise RuntimeError('depth > 60, LP bound stuck below 1? val=%s' % val)
    S = choose_cut(T, y, u, u_np, zero_edges)
    a, d = SIDES[S]
    return {'cut': S,
            'A': build(T + [a], u, u_np, zero_edges, stats, depth + 1),
            'D': build(T + [d], u, u_np, zero_edges, stats, depth + 1)}


def main():
    orbs = load_orbits(sys.argv[1])
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    names = sys.argv[3:] or sorted(orbs)
    for name in names:
        fn = os.path.join(out, name + '.json')
        if os.path.exists(fn):
            print(name, 'exists, skipped'); continue
        u = orbs[name]
        u_np = np.array([float(x) for x in u])
        zero_edges = [e for e in range(NE) if u[e] == 0]
        stats = {'nodes': 0, 'leaves': 0, 'depth': 0}
        t = time.time()
        tree = build([], u, u_np, zero_edges, stats)
        cert = {'orbit': name, 'u': [str(x) for x in u], 'tree': tree}
        with open(fn + '.tmp', 'w') as f:
            json.dump(cert, f)
        os.replace(fn + '.tmp', fn)
        print(f"{name}: zero edges {len(zero_edges)}, leaves {stats['leaves']}, depth {stats['depth']}, "
              f"time {time.time()-t:.1f}s", flush=True)


if __name__ == '__main__':
    main()
