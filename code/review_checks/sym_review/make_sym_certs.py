"""Generate rational branching certificates for the lower bound R*_det(v) >= U_det(v)
at every interior cell vertex v of the symmetric closed form (result 11, [S31]).

Logic of the certificate (checked independently by verify_sym_certs.py):
  * By result 2(b) the optimal deterministic rule can be taken to be a signed threshold
    rule; under symmetric parameters we may permute members so that its weights are
    non-increasing.  Such a rule F is self-dual and *regular*: moving one error from
    member j to member j-1 never turns a correct decision into an error.
  * A node of the tree fixes F on some complementary pairs (closed under regularity and
    self-duality).  A leaf carries a rational distribution pi in P(q,b); every rule F
    consistent with the node satisfies
        E_pi F >= sum_{x fixed to 1} pi_x + sum_{undecided pairs} min(pi_x, pi_{1-x}),
    and the leaf asserts that the right-hand side is >= U_det(v).
  * A node whose assignment is contradictory contains no regular self-dual rule.
Floats are used only to choose branches and candidate distributions; every leaf
distribution is rebuilt exactly with Fractions and checked before it is stored.

Usage: python3 -I make_sym_certs.py K [seconds]     (resumable; < 5 min per run)
Output: sym_certs/K{K}.json
"""
import os, sys, json, time, itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
_argv = sys.argv; sys.argv = [sys.argv[0], 'none']
import check_s31 as C            # closed-form terms and exact vertex enumeration
sys.argv = _argv

K = int(sys.argv[1]); TLIM = float(sys.argv[2]) if len(sys.argv) > 2 else 270
OUT = os.path.join(HERE, 'sym_certs'); os.makedirs(OUT, exist_ok=True)
PATH = os.path.join(OUT, f'K{K}.json')
T0 = time.time()

N = 1 << K; ALL = N - 1
E = list(itertools.combinations(range(K), 2))
REPS = [x for x in range(N) if x < (ALL ^ x)]
SUCC = [[x ^ (1 << j) ^ (1 << (j - 1)) for j in range(1, K) if (x >> j) & 1 and not (x >> (j - 1)) & 1] for x in range(N)]
PRED = [[] for _ in range(N)]
for x in range(N):
    for y in SUCC[x]:
        PRED[y].append(x)
BIT = np.array([[(x >> i) & 1 for i in range(K)] for x in range(N)])
PAIRM = np.array([[(x >> i) & (x >> j) & 1 for (i, j) in E] for x in range(N)])


def assign(val, x, v):
    """Set F(x)=v with self-duality and regularity propagation; return False on conflict."""
    stack = [(x, v)]
    while stack:
        y, w = stack.pop()
        if val[y] == w:
            continue
        if val[y] == 1 - w:
            return False
        val[y] = w
        stack.append((ALL ^ y, 1 - w))
        stack.extend((z, 1) for z in SUCC[y]) if w == 1 else stack.extend((z, 0) for z in PRED[y])
    return True


def node_lp(val, q, b):
    und = [r for r in REPS if val[r] < 0]
    nv = N + len(und)
    c = np.zeros(nv)
    for x in range(N):
        if val[x] == 1:
            c[x] = -1.0
    c[N:] = -1.0
    Aeq = np.vstack([np.ones(N), BIT.T]); Aeq = np.hstack([Aeq, np.zeros((Aeq.shape[0], len(und)))])
    beq = np.r_[1.0, np.full(K, q)]
    Aub = lil_matrix((len(E) + 2 * len(und), nv)); bub = np.zeros(len(E) + 2 * len(und))
    Aub[:len(E), :N] = PAIRM.T; bub[:len(E)] = b
    for k, r in enumerate(und):
        Aub[len(E) + 2 * k, N + k] = 1; Aub[len(E) + 2 * k, r] = -1
        Aub[len(E) + 2 * k + 1, N + k] = 1; Aub[len(E) + 2 * k + 1, ALL ^ r] = -1
    res = linprog(c, A_ub=csr_matrix(Aub), b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
    return -res.fun, res.x[:N], und


def exact_value(val, und, pi):
    v = sum((pi.get(x, Fr(0)) for x in range(N) if val[x] == 1), Fr(0))
    return v + sum((min(pi.get(r, Fr(0)), pi.get(ALL ^ r, Fr(0))) for r in und), Fr(0))


def exact_feasible(pi, q, b):
    if any(p < 0 for p in pi.values()) or sum(pi.values()) != 1:
        return False
    for i in range(K):
        if sum((p for x, p in pi.items() if (x >> i) & 1), Fr(0)) != q:
            return False
    for (i, j) in E:
        if sum((p for x, p in pi.items() if (x >> i) & (x >> j) & 1), Fr(0)) > b:
            return False
    return True


def solve_exact(rows, rhs, n, guess):
    """Gauss-Jordan over Fractions; free variables take rationalised guesses."""
    A = [r[:] + [h] for r, h in zip(rows, rhs)]
    piv = []; r0 = 0
    for col in range(n):
        p = next((i for i in range(r0, len(A)) if A[i][col] != 0), None)
        if p is None:
            continue
        A[r0], A[p] = A[p], A[r0]
        inv = 1 / A[r0][col]; A[r0] = [a * inv for a in A[r0]]
        for i in range(len(A)):
            if i != r0 and A[i][col] != 0:
                f = A[i][col]; A[i] = [a - f * c for a, c in zip(A[i], A[r0])]
        piv.append(col); r0 += 1
    if any(all(a == 0 for a in row[:n]) and row[n] != 0 for row in A):
        return None
    sol = [None] * n
    free = [c for c in range(n) if c not in piv]
    for c in free:
        sol[c] = guess[c]
    for i, col in enumerate(piv):
        sol[col] = A[i][n] - sum((A[i][c] * sol[c] for c in free), Fr(0))
    return sol


def rationalise(val, und, pif, q, b, U, tight_target):
    S = [x for x in range(N) if pif[x] > 1e-10]
    idx = {x: k for k, x in enumerate(S)}; n = len(S)
    rows, rhs = [], []
    rows.append([Fr(1)] * n); rhs.append(Fr(1))
    for i in range(K):
        rows.append([Fr((x >> i) & 1) for x in S]); rhs.append(q)
    for (i, j) in E:
        mom = sum(pif[x] for x in S if (x >> i) & (x >> j) & 1)
        if mom > float(b) - 1e-10:
            rows.append([Fr((x >> i) & (x >> j) & 1) for x in S]); rhs.append(b)
    obj = [Fr(0)] * n
    for x in S:
        if val[x] == 1:
            obj[idx[x]] += 1
    for r in und:
        a, c = r, ALL ^ r
        if a in idx and c in idx and abs(pif[a] - pif[c]) < 1e-10:
            row = [Fr(0)] * n; row[idx[a]] = Fr(1); row[idx[c]] = Fr(-1); rows.append(row); rhs.append(Fr(0))
            obj[idx[a]] += 1
        else:
            m = a if pif[a] <= pif[c] else c
            if m in idx:
                obj[idx[m]] += 1
    if tight_target:
        rows.append(obj); rhs.append(U)
    for den in (10 ** 6, 10 ** 9, 10 ** 12):
        guess = [Fr(float(pif[x])).limit_denominator(den) for x in S]
        sol = solve_exact(rows, rhs, n, guess)
        if sol is None:
            return None
        pi = {x: s for x, s in zip(S, sol) if s != 0}
        if exact_feasible(pi, q, b) and exact_value(val, und, pi) >= U:
            return pi
    return None


def build(val, q, b, U, stats, deadline):
    if time.time() > deadline:
        raise TimeoutError
    stats['nodes'] += 1
    lpv, pif, und = node_lp(val, q, b)
    if lpv >= float(U) - 1e-9:
        pi = rationalise(val, und, pif, q, b, U, tight_target=(lpv < float(U) + 1e-9))
        if pi is not None:
            stats['leaves'] += 1
            return {'pi': {str(x): str(p) for x, p in sorted(pi.items())}}
        if lpv > float(U) + 1e-7:
            raise RuntimeError('rationalisation failed with slack')
    if not und:
        raise RuntimeError(f'leaf below U: {lpv} < {float(U)}')
    r = max(und, key=lambda r: min(pif[r], pif[ALL ^ r]))
    node = {'pair': r, 'children': {}}
    for v in (1, 0):
        child = val[:]
        if not assign(child, r, v):
            stats['infeasible'] += 1
            node['children'][str(v)] = {'infeasible': True}
        else:
            node['children'][str(v)] = build(child, q, b, U, stats, deadline)
    return node


state = json.load(open(PATH)) if os.path.exists(PATH) else {'K': K, 'vertices': {}}
Tf = C.full_family(K); Tr = C.reduced_family(K)
verts = [(qq, bb, u) for qq, bb, u, tight, onb in C.cell_vertices(K, Tf) if u > 0 and qq < Fr(1, 2) and not onb]
print(f'K={K}: {len(verts)} interior vertices with U>0')
for qq, bb, u in verts:
    key = f'{qq},{bb}'
    if key in state['vertices']:
        continue
    stats = {'nodes': 0, 'leaves': 0, 'infeasible': 0}; t1 = time.time()
    try:
        tree = build([-1] * N, qq, bb, u, stats, T0 + TLIM)
    except TimeoutError:
        print(f'  time limit reached during vertex {key}; rerun to continue'); break
    state['vertices'][key] = {'q': str(qq), 'b': str(bb), 'U': str(u), 'stats': stats, 'tree': tree}
    json.dump(state, open(PATH, 'w'))
    print(f'  vertex ({key}) U={u}: {stats} [{time.time()-t1:.1f}s]', flush=True)
done = sum(1 for qq, bb, u in verts if f'{qq},{bb}' in state['vertices'])
print(f'K={K}: {done}/{len(verts)} vertices certified  [{time.time()-T0:.0f}s]')
