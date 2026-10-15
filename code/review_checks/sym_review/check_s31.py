"""Independent review of the symmetric-parameter closed form (uploaded report, archived as [S31]).

Usage: python3 -I check_s31.py PHASE K [seconds]
  PHASE = ident   : symbolic (sympy) check of the elimination identities for general s, t
          rand    : randomized closed form vs exact LP (random + boundary points)
          full    : reduced family min == full (uneliminated) family min on a rational grid
          det     : U_det vs exact deterministic MILP at (i) all cell vertices, enumerated
                    independently from the FULL family, (ii) boundary points, (iii) random points
          remark  : K=6 counterexample numbers (exact LP per rule)
Every run is < 5 minutes.
"""
import sys, time, itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

PHASE = sys.argv[1]
K = int(sys.argv[2]) if len(sys.argv) > 2 else 0
TLIM = float(sys.argv[3]) if len(sys.argv) > 3 else 240


# ---------- terms written directly from the report's formulas (a + B q + C b) ----------
def odd(s, t):   # B_{2s+1,t}
    d = (s + 1 - t) * (s - t)
    return (Fr(t * (t + 1), d), Fr(-2 * t * (2 * s + 1), d), Fr(2 * s * (2 * s + 1), d))


def even(s, t):  # B^even_{2s,t}
    d = (s - t) * (s - 1 - t)
    return (Fr(t * (t + 1), d), Fr(-2 * t * 2 * s, d), Fr(2 * s * (2 * s - 1), d))


def anch(s, t):  # B^anch_{2s,t}
    d = (s + 1 - t) * (s - t)
    return (Fr(t * (t + 1), d), Fr(2 - 2 * t * (2 * s + 1), d), Fr(2 * (2 * s - 1) * (s + 1), d))


def full_family(K):
    T = [(Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(3))]
    for s in range(1, (K - 1) // 2 + 1):
        for t in range(0, s):
            T.append(odd(s, t))
    for s in range(2, K // 2 + 1):
        for t in range(0, s - 1):
            if s - 1 - t > 0:
                T.append(even(s, t))
            T.append(anch(s, t))
    return sorted(set(T))


def reduced_family(K):
    T = [(Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(3))]
    if K % 2:
        m = (K - 1) // 2
        T += [odd(m, t) for t in range(1, m)]
        T += [odd(s, t) for s in range(2, m) for t in range(1, min(s - 1, (s - 1) // 3) + 1)]
    else:
        m = K // 2
        for t in range(1, m - 1):
            T += [odd(m - 1, t), even(m, t), anch(m, t)]
        T += [odd(s, t) for s in range(2, m - 1) for t in range(1, min(s - 1, (s - 1) // 3) + 1)]
    return sorted(set(T))


def U(T, q, b):
    return min(a + B * q + C * b for a, B, C in T)


def bmin(K, q):
    return max(Fr(2 * r * K, 1) * q - r * (r + 1) for r in range(K // 2 + 1)) / (K * (K - 1))


def feasible(K, q, b):
    return 0 <= q <= Fr(1, 2) and bmin(K, q) <= b <= q


# ---------- exact (float) solvers over all 2^K patterns ----------
def setup(K):
    E = list(itertools.combinations(range(K), 2))
    Z = np.array([[(s >> i) & 1 for i in range(K)] for s in range(1 << K)], float)
    P = np.array([[z[i] * z[j] for (i, j) in E] for z in Z])
    return E, Z, P


def rand_lp(K, q, b, E, Z, P):
    """randomized minimax = min L_Q s.t. Q(x)>=0, Q(x)+Q(1-x)>=1, gamma>=0 (dual LP)."""
    M = len(E); n = 1 + K + M; ALL = (1 << K) - 1
    c = np.r_[1.0, np.full(K, q), np.full(M, b)]
    A = []; lo = []
    for s in range(1 << K):
        A.append(np.r_[1.0, Z[s], P[s]]); lo.append(0.0)
        if s < (ALL ^ s):
            A.append(np.r_[2.0, Z[s] + Z[ALL ^ s], P[s] + P[ALL ^ s]]); lo.append(1.0)
    A = np.array(A); lo = np.array(lo)
    bnd = [(None, None)] * (1 + K) + [(0, None)] * M
    r = linprog(c, A_ub=-A, b_ub=-lo, bounds=bnd, method='highs')
    return r.fun


def det_milp(K, q, b, E, Z, P, sort_w=True, tl=120):
    M = len(E); ALL = (1 << K) - 1
    reps = [s for s in range(1 << K) if not (s >> (K - 1)) & 1]
    t = len(reps); n = 1 + K + M + t
    c = np.r_[1.0, np.full(K, q), np.full(M, b), np.zeros(t)]
    A = []; lo = []
    for j, r in enumerate(reps):
        row = np.zeros(n); row[0] = 1; row[1:1 + K] = Z[r]; row[1 + K:1 + K + M] = P[r]; row[1 + K + M + j] = -1
        A.append(row); lo.append(0.0)
        cc = ALL ^ r
        row = np.zeros(n); row[0] = 1; row[1:1 + K] = Z[cc]; row[1 + K:1 + K + M] = P[cc]; row[1 + K + M + j] = 1
        A.append(row); lo.append(1.0)
    if sort_w:  # w_i >= w_{i+1}, w_i = beta_i + 1/2 sum_j gamma_ij  (valid symmetry breaking)
        for i in range(K - 1):
            row = np.zeros(n); row[1 + i] = 1; row[2 + i] = -1
            for k, (u, v) in enumerate(E):
                if i in (u, v): row[1 + K + k] += .5
                if i + 1 in (u, v): row[1 + K + k] -= .5
            A.append(row); lo.append(0.0)
    lb = np.r_[np.full(1 + K, -np.inf), np.zeros(M + t)]; ub = np.r_[np.full(1 + K + M, np.inf), np.ones(t)]
    res = milp(c, integrality=np.r_[np.zeros(1 + K + M), np.ones(t)], bounds=Bounds(lb, ub),
               constraints=LinearConstraint(np.array(A), np.array(lo), np.inf),
               options={'time_limit': tl, 'mip_rel_gap': 1e-10})
    return res.fun, res.status


def cell_vertices(K, T):
    """Independent enumeration: all points of the polygon where the lower envelope of T has a vertex."""
    lines = [(Fr(1), Fr(0), Fr(0)), (Fr(1), Fr(0), Fr(1, 2)), (Fr(1), Fr(-1), Fr(0))]  # q=0, q=1/2, q-b=0
    for r in range(K // 2 + 1):  # 2rKq - K(K-1) b = r(r+1)
        lines.append((Fr(2 * r * K), Fr(-K * (K - 1)), Fr(r * (r + 1))))
    for (a1, B1, C1), (a2, B2, C2) in itertools.combinations(T, 2):
        lines.append((B1 - B2, C1 - C2, a2 - a1))
    pts = set()
    for (A1, B1, D1), (A2, B2, D2) in itertools.combinations(lines, 2):
        det = A1 * B2 - A2 * B1
        if det == 0:
            continue
        q = (D1 * B2 - D2 * B1) / det; b = (A1 * D2 - A2 * D1) / det
        if feasible(K, q, b):
            pts.add((q, b))
    out = []
    for q, b in pts:
        u = U(T, q, b)
        tight = sum(1 for a, B, C in T if a + B * q + C * b == u)
        onb = (q == 0) + (q == Fr(1, 2)) + (b == q) + (b == bmin(K, q))
        corner = any(q == Fr(r, K) for r in range(K // 2 + 1)) and b == bmin(K, q)
        if tight >= 3 or (tight >= 2 and onb) or (onb and corner) or (q in (0, Fr(1, 2)) and b in (q, bmin(K, q))):
            out.append((q, b, u, tight, bool(onb)))
    return sorted(out)


T0 = time.time()
if PHASE == 'ident':
    import sympy as sp
    s, t, q, b = sp.symbols('s t q b')
    Bo = lambda S, T_: (T_ * (T_ + 1) - 2 * T_ * (2 * S + 1) * q + 2 * S * (2 * S + 1) * b) / ((S + 1 - T_) * (S - T_))
    Be = (t * (t + 1) - 2 * t * 2 * s * q + 2 * s * (2 * s - 1) * b) / ((s - t) * (s - 1 - t))
    Ba = (t * (t + 1) + (2 - 2 * t * (2 * s + 1)) * q + 2 * (2 * s - 1) * (s + 1) * b) / ((s + 1 - t) * (s - t))
    e1 = sp.simplify(Ba - Bo(s, t) - 2 * (q - b) / ((s + 1 - t) * (s - t)))
    c1 = 1 - t / (2 * s + 1) + 2 * s / ((2 * s + 1) * (s - 1 - t)); c2 = t / (2 * s + 1)
    e2 = sp.simplify(Be - c1 * Bo(s, t) - c2 * Bo(s, t + 1))
    L1 = (2 * s - 2 * t - 3) * (s - t + 1) / ((2 * s + 1) * (s - t - 1)); L2 = 2 * t / (2 * s + 1); L3 = sp.Rational(2, 3) / ((s - t) * (s - t - 1))
    e3 = sp.simplify(Bo(s - 1, t) - L1 * Bo(s, t) - L2 * Bo(s, t + 1) - L3 * 3 * b)
    e4 = sp.simplify(L1 + L2 + L3 - 1 - 2 * (3 * t + 1 - s) / (3 * (2 * s + 1) * (s - t) * (s - t - 1)))
    print('identity residuals (should all be 0):', e1, e2, e3, e4)
    # sign conditions of L1 in the stated range 1 <= t <= s-2
    print('L1 numerator factor (2s-2t-3) >= 1 for t <= s-2:', all((2 * S - 2 * T_ - 3) >= 1 for S in range(3, 60) for T_ in range(1, S - 1)))
    # randomized closed form: D formula
    for KK in range(3, 30):
        m = KK // 2
        for tt in range(0, (KK + 1) // 2 - 1):
            D = (m - tt) * (m - tt - 1) + (KK - m - tt) * (KK - m - tt - 1)
            assert D == (2 * (m - tt) ** 2 if KK % 2 else 2 * (m - tt) * (m - tt - 1)) and D > 0
    print('D^rand formula and positivity: OK for K = 3..29')

elif PHASE == 'remark':
    E, Z, P = setup(6); q, b = 5 / 18, 1 / 18
    def worst(Fv):
        A_eq = np.vstack([np.ones(64), Z.T]); r = linprog(-Fv, A_ub=P.T, b_ub=np.full(15, b), A_eq=A_eq, b_eq=np.r_[1, np.full(6, q)], bounds=(0, None), method='highs'); return -r.fun
    for w in [(1, 0, 0, 0, 0, 0), (1, 1, 1, 0, 0, 0), (1, 1, 1, 1, 1, 0), (2, 1, 1, 1, 1, 1)]:
        w = np.array(w, float); Fv = (Z @ w > w.sum() / 2).astype(float)
        print('rule', tuple(int(x) for x in w), 'worst-case =', round(worst(Fv), 10))
    print('exact R*_det =', det_milp(6, q, b, E, Z, P)[0], ' claimed B^anch_{6,1} = 4/27 =', 4 / 27)

elif PHASE == 'rand':
    E, Z, P = setup(K); rng = np.random.default_rng(100 + K); m = K // 2
    def rand_cf(q, b):
        v = [q]
        for t in range(0, (K + 1) // 2 - 1):
            D = (m - t) * (m - t - 1) + (K - m - t) * (K - m - t - 1)
            v.append((t * (t + 1) - 2 * t * K * q + K * (K - 1) * b) / D)
        return min(v)
    n = 0; worst = 0
    while time.time() - T0 < TLIM and n < 400:
        q = Fr(int(rng.integers(1, 1000)), 2000)
        mode = n % 3
        if mode == 0: b = bmin(K, q)                     # lower boundary
        elif mode == 1: b = bmin(K, q) + (q - bmin(K, q)) * Fr(int(rng.integers(0, 1000)), 1000)
        else: b = bmin(K, q) + (q - bmin(K, q)) * Fr(int(rng.integers(0, 60)), 1000)  # near lower boundary
        ex = rand_lp(K, float(q), float(b), E, Z, P); cf = rand_cf(float(q), float(b))
        worst = max(worst, abs(ex - cf)); n += 1
    print(f'K={K} randomized closed form vs exact LP: {n} points, max |diff| = {worst:.2e}  [{time.time()-T0:.0f}s]')

elif PHASE == 'full':
    Tf = full_family(K); Tr = reduced_family(K); worst = Fr(0); n = 0; bad = 0
    for i in range(0, 201):
        q = Fr(i, 400)
        lo = bmin(K, q)
        for j in range(0, 41):
            b = lo + (q - lo) * Fr(j, 40)
            uf, ur = U(Tf, q, b), U(Tr, q, b); n += 1
            if ur != uf: bad += 1; worst = max(worst, ur - uf)
    print(f'K={K}: reduced family ({len(Tr)} terms) vs full family ({len(Tf)} terms) on {n} rational grid points: mismatches = {bad}, max diff = {worst}')

elif PHASE == 'det':
    E, Z, P = setup(K); Tf = full_family(K); Tr = reduced_family(K)
    V = cell_vertices(K, Tf)
    nontriv = [v for v in V if v[2] > 0 and v[0] < Fr(1, 2)]
    print(f'K={K}: {len(V)} cell vertices of the full-family envelope ({len(nontriv)} with U>0 and q<1/2)')
    worst = 0; nchk = 0
    for q, b, u, tight, onb in nontriv:
        assert U(Tr, q, b) == u
        val, st = det_milp(K, float(q), float(b), E, Z, P)
        d = val - float(u); worst = max(worst, abs(d)); nchk += 1
        flag = '' if abs(d) < 1e-8 else '  <-- MISMATCH'
        print(f'   ({q}, {b})  U={u}  MILP={val:.9f}  diff={d:.1e}  tight={tight} boundary={onb}{flag}')
        if time.time() - T0 > TLIM: print('   time limit reached'); break
    rng = np.random.default_rng(K); nr = 0; wr = 0
    while time.time() - T0 < TLIM and nr < 40:
        q = Fr(int(rng.integers(1, 500)), 1000); lo = bmin(K, q)
        b = lo + (q - lo) * Fr(int(rng.integers(0, 1000)), 1000) * (Fr(1, 4) if nr % 2 else 1)
        val, st = det_milp(K, float(q), float(b), E, Z, P)
        wr = max(wr, abs(val - float(U(Tr, q, b)))); nr += 1
    print(f'K={K}: vertices checked {nchk}/{len(nontriv)}, max |diff| = {worst:.1e}; random points {nr}, max |diff| = {wr:.1e}  [{time.time()-T0:.0f}s]')
