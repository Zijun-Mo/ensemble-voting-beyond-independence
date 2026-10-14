"""Brute-force count of self-dual monotone Boolean functions on K variables (K<=5) and check that all are
positive-weight threshold functions (i.e. our LTF enumeration is complete)."""
import itertools, numpy as np
from core import patterns, flip_index, reps, enumerate_selfdual_ltfs
for K in (3, 4, 5):
    Z = patterns(K).astype(int); fl = flip_index(K); R = reps(K); n = 2 ** K
    # covering pairs: v -> v + e_i
    edges = [(v, v | (1 << i)) for v in range(n) for i in range(K) if not v & (1 << i)]
    cnt = 0; tabs = set()
    for bits in itertools.product((0, 1), repeat=len(R)):
        g = np.zeros(n, int)
        for b, v in zip(bits, R):
            g[v] = b; g[fl[v]] = 1 - b
        if all(g[a] <= g[b] for a, b in edges):
            cnt += 1; tabs.add(tuple(g))
    G, _ = enumerate_selfdual_ltfs(K, n_random=20000)
    ltf = set(map(tuple, G))
    print(f"K={K}: self-dual monotone functions = {cnt}; LTFs found = {len(ltf)}; all monotone self-dual are LTF: {tabs == ltf}")
