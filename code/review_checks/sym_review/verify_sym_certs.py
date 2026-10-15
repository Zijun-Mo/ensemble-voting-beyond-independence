"""Independent exact verifier for the branching certificates in sym_certs/K*.json.

Standard library only (fractions, itertools, json).  Nothing is imported from the
generator.  For each K in 3..11 it

 1. rebuilds the symmetric closed form U_det from the formulas of [S31] (dictator, 3b,
    B_{2s+1,t}, B^even_{2s,t}, B^anch_{2s,t}; full uneliminated family),
 2. enumerates exactly every vertex of the cell decomposition of U_det over the
    feasible polygon Theta^sym_K and checks that each vertex that is interior to the
    polygon and has U_det > 0 carries a certificate,
 3. walks every certificate tree: a branch fixes F on one complementary pair; the
    assignment is closed under self-duality F(1-x) = 1 - F(x) and regularity
    (moving an error from member j to member j-1 cannot turn a correct decision into
    an error); an 'infeasible' leaf must yield a contradiction; a distribution leaf
    must be exactly feasible (marginals = q, pair moments <= b) and satisfy
        sum_{F(x)=1 fixed} pi_x + sum_{undecided pairs} min(pi_x, pi_{1-x}) >= U_det(v).

Together with result 2(b) (threshold sufficiency; sorting the weights makes the rule
regular), concavity of each rule's worst-case risk in (q, b), and the analytic
boundary theorem of [S31], this proves R*_det = U_det on the whole polygon for K <= 11.
Usage: python3 -I verify_sym_certs.py
"""
import json, os, itertools, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def terms(K):
    T = {(Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(3))}
    for s in range(1, (K - 1) // 2 + 1):          # odd committees of size 2s+1 <= K
        for t in range(s):
            d = (s + 1 - t) * (s - t)
            T.add((Fr(t * t + t, d), Fr(-2 * t * (2 * s + 1), d), Fr(2 * s * (2 * s + 1), d)))
    for s in range(2, K // 2 + 1):                 # even layer / tie-break, size 2s <= K
        for t in range(s - 1):
            if s - 1 - t > 0:
                d = (s - t) * (s - 1 - t)
                T.add((Fr(t * t + t, d), Fr(-4 * s * t, d), Fr(2 * s * (2 * s - 1), d)))
            d = (s + 1 - t) * (s - t)
            T.add((Fr(t * t + t, d), Fr(2 - 2 * t * (2 * s + 1), d), Fr(2 * (2 * s - 1) * (s + 1), d)))
    return sorted(T)


def lower(K, q):
    return max((2 * h * K * q - h * (h + 1)) / Fr(K * (K - 1)) for h in range(K // 2 + 1))


def interior_vertices(K):
    T = terms(K)
    U = lambda q, b: min(a + B * q + c * b for a, B, c in T)
    bnd = [(Fr(1), Fr(0), Fr(0)), (Fr(1), Fr(0), Fr(1, 2)), (Fr(1), Fr(-1), Fr(0))]
    bnd += [(Fr(2 * h * K), Fr(-K * (K - 1)), Fr(h * (h + 1))) for h in range(K // 2 + 1)]
    diff = [(B1 - B2, c1 - c2, a2 - a1) for (a1, B1, c1), (a2, B2, c2) in itertools.combinations(T, 2)
            if (B1, c1) != (B2, c2)]
    pts = set()
    for (A1, B1, D1), (A2, B2, D2) in itertools.combinations(bnd + diff, 2):
        det = A1 * B2 - A2 * B1
        if det:
            q = (D1 * B2 - D2 * B1) / det; b = (A1 * D2 - A2 * D1) / det
            if 0 <= q <= Fr(1, 2) and lower(K, q) <= b <= q:
                pts.add((q, b))
    out = []
    for q, b in pts:
        u = U(q, b)
        on_boundary = q in (0, Fr(1, 2)) or b == q or b == lower(K, q)
        ntight = sum(1 for a, B, c in T if a + B * q + c * b == u)
        # an interior point is a cell vertex iff at least three affine pieces are active
        if not on_boundary and ntight >= 3 and u > 0:
            out.append((q, b, u))
    return sorted(out)


def check_tree(K, q, b, U, tree):
    N = 1 << K; ALL = N - 1
    pairs = list(itertools.combinations(range(K), 2))

    def close(val):
        """Fixpoint of self-duality and regularity; returns None on contradiction."""
        changed = True
        while changed:
            changed = False
            for x in range(N):
                v = val.get(x)
                if v is None:
                    continue
                c = ALL ^ x
                if val.get(c, 1 - v) != 1 - v:
                    return None
                if c not in val:
                    val[c] = 1 - v; changed = True
                for j in range(1, K):
                    if (x >> j) & 1 and not (x >> (j - 1)) & 1:
                        y = x ^ (1 << j) ^ (1 << (j - 1))     # error moved from j to j-1
                        if v == 1 and val.get(y) != 1:
                            if val.get(y) == 0:
                                return None
                            val[y] = 1; changed = True
                        if v == 0 and val.get(y) == 0:
                            pass
                    if (x >> (j - 1)) & 1 and not (x >> j) & 1 and v == 0:
                        y = x ^ (1 << j) ^ (1 << (j - 1))     # x is obtained from y by a left move
                        if val.get(y) != 0:
                            if val.get(y) == 1:
                                return None
                            val[y] = 0; changed = True
        return val

    counts = {'leaves': 0, 'infeasible': 0}

    def walk(node, val):
        if 'infeasible' in node:
            assert close(dict(val)) is None, 'branch marked infeasible is consistent'
            counts['infeasible'] += 1
            return
        if 'pi' in node:
            pi = {int(x): Fr(p) for x, p in node['pi'].items()}
            assert all(p >= 0 for p in pi.values()) and sum(pi.values()) == 1
            for i in range(K):
                assert sum(p for x, p in pi.items() if (x >> i) & 1) == q
            for i, j in pairs:
                assert sum(p for x, p in pi.items() if (x >> i) & (x >> j) & 1) <= b
            bound = sum((p for x, p in pi.items() if val.get(x) == 1), Fr(0))
            bound += sum((min(pi.get(x, Fr(0)), pi.get(ALL ^ x, Fr(0)))
                          for x in range(N) if x < (ALL ^ x) and x not in val), Fr(0))
            assert bound >= U, (bound, U)
            counts['leaves'] += 1
            return
        r = node['pair']
        for v in (1, 0):
            child = dict(val); child[r] = v
            closed = close(child)
            sub = node['children'][str(v)]
            if closed is None:
                assert 'infeasible' in sub, 'contradictory branch not marked'
                counts['infeasible'] += 1
            else:
                walk(sub, closed)

    root = close({})
    walk(tree, root)
    return counts


T0 = time.time(); total = 0
for K in range(3, 12):
    verts = interior_vertices(K)
    path = os.path.join(HERE, 'sym_certs', f'K{K}.json')
    data = json.load(open(path))['vertices'] if os.path.exists(path) else {}
    for q, b, u in verts:
        key = f'{q},{b}'
        assert key in data, f'K={K}: vertex {key} has no certificate'
        rec = data[key]
        assert Fr(rec['U']) == u and Fr(rec['q']) == q and Fr(rec['b']) == b
        c = check_tree(K, q, b, u, rec['tree'])
        total += 1
        print(f'K={K:2d} vertex ({q}, {b})  U_det = {u}:  {c["leaves"]} distribution leaves, '
              f'{c["infeasible"]} contradictory branches  -> R*_det >= U_det certified')
    if not verts:
        print(f'K={K:2d}: no interior vertex with U_det > 0 (boundary theorem suffices)')
print(f'\n{total} vertices certified exactly in {time.time() - T0:.0f}s')
