#!/usr/bin/env python3
# Utilities for counting hypergraphs in subsets of vertices.

import math
import pdb
import sys

import numpy as np
import scipy.optimize
import scipy.sparse
import scipy.special

def log_factorial(n):
    """Natural log of factorial."""
    return math.lgamma(n+1)

def log_comb(n, k):
    """Computes log(n \choose k)."""
    return log_factorial(n) - (log_factorial(k) + log_factorial(n-k))

def log_sum_exp(a, b):
    """Computes log( exp(a) + exp(b) ), where a and b are vectors.

    Wrapper for logsumexp.
    """
    x = np.stack([a, b], axis=1)
    return scipy.special.logsumexp(x, b=[1,1], axis=1)

def log_diff_exp(b, a):
    """Computes log( exp(b) - exp(a) ), where a and b are vectors."""
    x = np.stack([b, a], axis=1)
    return scipy.special.logsumexp(x, b=[1,-1], axis=1)

class LogHypergraphCounter:
    """Counts hypergraphs in subsets of vertices (in log-space).

    When we restrict to a subset of vertices, presumably we'll
    need to keep track of how many hypergraphs have a given
    number of hyperedges.

    Presumably there's not a simpler expression for this...
    """
    def __init__(self, n, k):
        """Constructor.
   
        n: number of vertices in the larger graph
        k: size of cliques (hyperedges) 
        """
        self.n = n
        self.k = k

    def count_hypergraphs_max_vertices(self):
        """Counts hypergraphs with _up to_ some number of vertices.

        Returns: a dict h, with key v, an int in the range k..n,
            representing the number of vertices.
            h[v] is a numpy array of length ${v \choose k}$.
            h[v][i] is the number of hypergraphs with up to
                $v$ vertices, and exactly $i$ hyperedges.
        """
        # to compute this, first compute number of hypergraphs with
        # exactly some number of vertices used
        exact_counts = self.count_hypergraphs_exact_vertices()
        # this will hold the vectors of counts, for each number of vertices
        h = dict()
        # then, loop through the number of vertices
        for i in range(self.k, self.n+1):
            # start with a count of 0 (or -inf, in log-space)
            num_cliques = scipy.special.comb(i, self.k, exact=True)
            h[i] = np.full([num_cliques + 1], -np.inf)
            # add in number of hypergraphs with up to that many vertices
            for j in range(self.k, i+1):
                n1 = exact_counts[j].shape[0]
                h[i][:n1] = log_sum_exp(h[i][:n1],
                    log_comb(self.n, j) + exact_counts[j])
            # lastly, count the empty hypergraph (1, or 0 in log-space)
            h[i][0] = 0.
        return h

    def count_hypergraphs_exact_vertices(self):
        """Counts hypergraphs with _exactly_ some number of vertices.

        The nice thing about these is that since they all "use" some
        number of vertices, they're all distinct.
        This is used by "count_hypergraphs_max_vertices".
        """
        exact_counts = dict()
        # loop through number of vertices
        for i in range(self.k, self.n+1):
            # start with count of hypergraphs (on all n vertices) with
            # _up to_ this many vertices
            num_cliques = scipy.special.comb(i, self.k, exact=True)
            # note that this is in log-space
            exact_counts[i] = np.array([log_comb(num_cliques, r)
                for r in range(num_cliques + 1)])
            # also, don't count case with zero hypergraphs (as that's not
            # specific to a particular vertex set). Technically, this should
            # be -inf, but I don't think it's used anyway.
            exact_counts[i][0] = 0
            # then, subtract off hypergraphs with fewer vertices (if any)
            for j in range(self.k, i):
                n1 = exact_counts[j].shape[0]
                exact_counts[i][:n1] = log_diff_exp(exact_counts[i][:n1],
                    log_comb(i, j) + exact_counts[j])
        return exact_counts

# XXX a quick test
if __name__ == '__main__':
    # pdb.set_trace()

    # a small example
    hc = LogHypergraphCounter(4, 2)

    print('exact-number-of-vertex counts:')
    for (v, h) in hc.count_hypergraphs_exact_vertices().items():
        h = h.astype(int).tolist()
        print(f'{v}: {h}\n')

    print('up-to-some-number-of-vertex counts:')
    for (v, h) in hc.count_hypergraphs_max_vertices().items():
        h = h.astype(int).tolist()
        print(f'{v}: {h}\n')

