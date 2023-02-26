#!/usr/bin/env python3
# Utilities for counting hypergraphs in subsets of vertices.

import math
import pdb
import sys

import numpy as np
import scipy.optimize
import scipy.sparse
import scipy.special

class HypergraphCounter:
    """Counts hypergraphs in subsets of vertices.

    When we restrict to a subset of vertices, presumably we'll
    need to keep track of how many hypergraphs have a given
    number of hyperedges.
    """
    def __init__(self, n, k):
        """Constructor.
   
        n: number of vertices in the larger graph
        k: size of cliques (hyperedges) 
        """
        self.n = n
        self.k = k



    def count_hypergraphs_max_vertices(self):
        """

        Returns: a hash h, with key v, an int in the range k..n,
            representing the number of vertices.
            h[v] is a numpy array of length ${v \choose k}$.
            h[v][i] is the number of hypergraphs with up to
                $v$ vertices, and exactly $i$ hyperedges.

        """
        pass


    def count_hypergraphs_exact_vertices(self):
        """Counts hypergraphs with _exactly_ some num. of vertices.

        The nice thing about these, is that they're non-overlapping.
        This is used by "count_hypergraphs_max_vertices".
        """
        pass

def count_bits(x):
    """Counts number of bits set in a numpy vector."""
    # we assume numbers are "somewhat small" (if they were
    # large, they'd take a while to loop through anyway)
    num_bits_set = np.zeros(x.shape[0])
    for i in range(30):
        mask = np.array(2**i)
        num_bits_set += ((x & mask) > 0) + 0
    return num_bits_set

class SlowHypergraphCounter:
    """Slower reference implementation of HypergraphCounter.

    """
    def __init__(self, n, k):
        """ Constructor."""
        # restrict problem size
        if scipy.special.comb(n, k, exact=True) >= 40:
            raise ValueError(
                'n choose k would use a lot of memory; exiting')
        self.n = n
        self.k = k
        # numbering for the cliques
        cliques = [frozenset(s)
            for s in itertools.combinations(range(n), k)]
        self.clique_index = dict(zip(cliques, range(len(cliques))))
        # the sets of cliques
        self.S = np.arange(self.num_functions)
        # the sets of cliques which are "zeroable", and don't overlap S
        self.Z = np.zeros(self.num_functions, dtype='int')

    def count_zero_set(self, cliques_to_count):
        """Adds everywhere that some cliques could be zeroed.

        cliques_to_count: the cliques to count, as an int
            representing a bitset
        Side effects: fills in bits of Z which correspond to cliques
            which aren't in this set (in S), but could be zeroed
        """
        # find sets where those cliques would be missed
        i = (self.S & cliques_to_count) == 0
        # add those to Z
        self.Z[i] = self.Z[i] | cliques_to_count

    def zero_vertices(self):
        """Finds sets of cliques unaffected by zeroing out vertices.

        This finds the sets of cliques which wouldn't be 'hit'
        by zeroing out edges incident to some vertex.
        """
        # loop through the vertices
        for v in range(self.n):
            # find cliques including that vertex
            mask = 0
            vertex_set = frozenset([v])
            # loop through the cliques
            for (clique, i) in self.clique_index.items():
                # if the edge is in this clique, set the appropriate bit
                if clique > vertex_set:
                    mask |= 2 ** i
            # zero out that set of cliques
            self.count_zero_set(mask)

    def count (self):
        """


