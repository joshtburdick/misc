#!/usr/bin/env python3
# Attempt based on zeroing out edges (or vertices).

import pdb

import itertools
import numpy as np
import scipy.special
import scipy.stats

import lp_helper

def count_bits(x):
    """Counts number of bits set in a numpy vector."""
    # we assume numbers are "somewhat small" (if they were
    # large, they'd take a while to loop through anyway)
    num_bits_set = np.zeros(x.shape[0])
    for i in range(30):
        mask = np.array(2**i)
        num_bits_set += ((x & mask) > 0) + 0
    return num_bits_set

class LpBruteForce1:
    """Attempt at bound by finding zeroable edges/vertices, using brute force.

    """

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        This will have one variable for each possible number of cliques
        (or "level"). 

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        # the number of possible cliques
        self.num_cliques = scipy.special.comb(n, k, exact=True)
        # the number of functions
        self.num_functions = 2 ** self.num_cliques
        # wrapper for LP solver, with one variable per "level"
        self.lp = lp_helper.LP_Helper(range(self.num_cliques+1))
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
        """
        # find sets where those cliques would be missed
        i = (self.S & cliques_to_count) == 0
        # add those to Z
        self.Z[i] = self.Z[i] | cliques_to_count

    def cliques_hit_by_edges(self, edge):
        """Gets the set of cliques hit by an edge.

        edge: the edge (as a pair)
        Returns: the cliques hit by that edge, as an int
            with bits set according to self.clique_index.
        """
        edge1 = frozenset(edge)
        # loop through the cliques
        mask = 0
        # loop through the cliques
        for (clique, i) in self.clique_index.items():
            # if the edge is in this clique, set the appropriate bit
            if clique > edge1:
                mask |= 2 ** i
        return mask

    def zero_edges(self):
        # loop through the edges
        for edge in itertools.combinations(range(self.n), 2):
            # record which sets of cliques are "missed"
            mask = self.cliques_hit_by_edges(edge)
            self.count_zero_set(mask)

    def zero_vertices(self):
        """Presumably gives a worse (but easier to analyze) bound."""

        pass

    def add_upper_bound_constraints(self):
        """Adds constraints based on upper bounds from zeroing."""
        # count cliques in each set
        set_size = count_bits(self.S)
        # count "zeroable" cliques "above" each set
        zeroable_size = count_bits(self.Z)
        # loop through the "levels" (number of cliques in C)
        for i in range(self.num_cliques+1):
            # find sets with this many cliques
            j = (set_size == i)
            # these should just be binomial coefficients
            assert(sum(j)
                == scipy.special.comb(self.num_cliques, i, exact=True) )
            # get average number of zeroable sets "above" these
            mean_zeroable = np.mean(zeroable_size[j])
            print(str(i) + '   ' + str(mean_zeroable))
            # add constraint that a function at this level...
            self.lp.add_constraint([(i, 1.)],
                '<',
                # ... has expected rank less than the number zeroable,
                # minus half the number of sets of this size
                mean_zeroable - np.sum(j) / 2.)

    def add_average_rank_constraint(self):
        """Adds equality constraint on average rank of all of the functions."""
        self.lp.add_constraint(
            [(i, scipy.special.comb(self.num_cliques, i, exact=False) / self.num_functions)
                for i in range(self.num_cliques+1)],
            '=',
            self.num_functions / 2.)

    def get_all_bounds(self):
        """Gets bounds for all the sets.

        This solves the linear program, including all constraints which
            have been included.
        Returns: a hash, of the bounds for all of the variables.
        """
        # we're trying to see how low the rank of
        # "finding all the cliques" can go
        bounds = self.lp.solve(self.num_cliques,
            bounds=(0, self.num_functions-1))
        return bounds

if __name__ == '__main__':
    bound = LpBruteForce1(4, 3)
    bound.add_average_rank_constraint()
    bound.zero_edges()
    bound.add_upper_bound_constraints()
    pdb.set_trace()
    b = bound.get_all_bounds()
    for i in range(bound.num_cliques+1):
        print('\t'.join([str(i), str(b[i])]))
    print()
    print('bound = ' + str(b[bound.num_cliques]))

