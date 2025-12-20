#!/usr/bin/env python3
# Checks what fraction of hypergraphs which cover all vertices
# have a given number of hyperedges.

import sys

sys.path.append('..')

from scipy.special import comb
import hypergraph_counter

def fraction_covering_hypergraphs(n, k):
    """Checks what fraction of hypergraphs which cover all vertices 
    have "a lot of hyperedges.

    In particular, it's convenient to define "a lot of hyperedges" as 
    at least comb(n-1, k)+1 hyperedges, since then, the hypergraph must
    cover all vertices.
    
    n, k: sizes of hypergraphs to consider
    Returns: fraction of hypergraphs which cover all vertices
    """
    counter = hypergraph_counter.HypergraphCounter(n, k)
    function_counts = counter.count_hypergraphs_exact_vertices()[n]
    total = function_counts.sum()
    threshold = comb(n-1, k, exact=True) + 1
    covering = function_counts[threshold:].sum()
    return {
        'n': n,
        'k': k,
        'total': total,
        'covering': covering,
        'threshold': threshold,
        'fraction': covering / total
    }



for n in range(5, 14):
    k = n // 2 + 1
    stats = fraction_covering_hypergraphs(n, k)
    print(f"n={n}, k={k}: {stats['fraction']:.3f}")
    

sys.exit(0)



if len(sys.argv) != 4:
    print("Usage: python3 rank_check.py <n_low> <n_high> <k>")
    sys.exit(1)

n_low = int(sys.argv[1])
n_high = int(sys.argv[2])
k = int(sys.argv[3])

for n in range(n_low, n_high+1):
    fraction = fraction_covering_hypergraphs(n, k)
    print(f"n={n}, k={k}: {fraction:.3f}")
