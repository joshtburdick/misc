#!/usr/bin/env python3
# Utilities for counting hypergraphs in subsets of vertices.

import math
import pdb
import sys

import numpy as np
import scipy.optimize
import scipy.sparse
import scipy.special


def log_comb(n, k):
    """Computes log(n \choose k).
    
    """
    return math.lgamma(n) - (math.lgamma(k) + math.lgamma(n-k)

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
        """Counts hyperghraphs with _up to_ some number of vertices.

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
            # start with a count of 0
            num_cliques = scipy.special.comb(i, self.k, exact=True)
            h[i] = np.full([num_cliques + 1], 0.)
            # add in number of hypergraphs with up to that many vertices
            for j in range(self.k, i+1):
                n1 = exact_counts[j].shape[0]
    # FIXME this should use numpy.logaddexp
                h[i][:n1] += scipy.special.comb(self.n, j, exact=True) * exact_counts[j]
            # lastly, count the empty hypergraph
            h[i][0] = 1
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
            exact_counts[i] = np.array([scipy.special.comb(num_cliques, r, exact=False)
                for r in range(num_cliques + 1)])
            # also, don't count case with zero hypergraphs (as that's not
            # specific to a particular vertex set)
            exact_counts[i][0] = 0
            # then, subtract off hypergraphs with fewer vertices (if any)
            for j in range(self.k, i):
                n1 = exact_counts[j].shape[0]
                exact_counts[i][:n1] -= scipy.special.comb(i, j, exact=False) * exact_counts[j]
        return exact_counts

# XXX a quick test
if __name__ == '__main__':
    # FIXME not yet implemented
    return
    # this is a toy example, but is small enough to check by hand
    # hc = HypergraphCounter(4, 2)
    # a slightoly larger example
    hc = HypergraphCounter(6, 3)

    print('exact-number-of-vertex counts:')
    for (v, h) in hc.count_hypergraphs_exact_vertices().items():
        h = h.astype(int).tolist()
        print(f'{v}: {h}\n')

    print('up-to-some-number-of-vertex counts:')
    for (v, h) in hc.count_hypergraphs_max_vertices().items():
        h = h.astype(int).tolist()
        print(f'{v}: {h}\n')

