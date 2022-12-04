#!/usr/bin/env python3
# Attempt at bounding the average rank of functions,
# relative to the number of hyperedges.

import pdb

from scipy.special import comb

import lp_helper

class RankBound3:
    """Bound, using the lattice of canonical hypergraphs.

    """

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        """
        # the number of possible cliques
        self.num_cliques = comb(n, k, exact=True)
        # the number of functions
        self.num_functions = 2 ** self.num_cliques
        # number of cliques which could be "hit" by feeding in a 0 to an edge
        self.max_cliques_zeroed = comb(n-2, k-2, exact=True)
        # Set up linear program. For each "level" (number of hyperedges),
        # there are three sets, with associated variables:
        #   A: sets which are hit by edge e
        #   B: sets which are missed by edge e
        #   C: A \union B
        # The variables track the expected rank of functions finding those
        # sets of hyperedges.
        variables = ([('A', i) for i in range(self.num_cliques+1)]
            + [('B', i) for i in range(self.num_cliques+1)]
            + [('C', i) for i in range(self.num_cliques+1)])
        # wrapper for LP solver
        self.lp = lp_helper.LP_Helper(variables)

    def add_mixture_equalities(self):
        """Add bound that |C| is a mixture of |A| and |B|."""
        for i in range(self.num_cliques+1):
            max_z = min(i, self.max_cliques_zeroed)
            # probability that feeding in a 0 misses all the cliques
            b = scipy.special.hyperg(self.num_cliques, i, max_z)(0)
            self.lp.add_constraint(
                [(('A', i), 1-b), (('B', i), b), (('C', i), -1.)],
                '=', 0.)
            # ??? can this be folded into add_zeroing_constraints() ?

    def add_zeroing_constraints(self):
        """Adds constraints from zeroing out one edge."""
        # loop through the "levels" (number of cliques in C)
        for i in range(self.num_cliques+1):
            # this is the maximum number of cliques which could
            # potentially be zeroed (at this level)
            max_z = min(i, self.max_cliques_zeroed)
            # the probability of some number of edges being zeroed
            z = scipy.special.hyperg(self.num_cliques, i, max_z)
            self.lp.add_constraint(
                # constraint that A (which is "hit") is a weighted sum of
                # what's left in B (none of which are "hit"), after zonking
                [('A', i), 1.)] + [
                    (('B', i-j), -z(j))
                    for j in range(max_z+1)],
                '>=',
                # the expected value of this is at least half the number
                # of functions being added in A
                0.5 * (1. - z(0)) * comb(self.num_cliques, i, exact=True)

    def add_average_rank_constraint(self):
        """Adds equality constraint on average rank of all of the functions."""
        self.lp.add_constraint(
            [(('C', i), comb(self.num_cliques, i, exact=True))
                for i in range(self.num_cliques+1)],
            '=',
            (2**self.num_cliques) / 2.)

    def get_all_bounds(self):
        """Gets bounds for all the sets.

        This solves the linear program, including all constraints which
            have been included.
        Returns: a hash, of the bounds for all of the variables.
        """
        # we're trying to see how low the rank of
        # "finding all the cliques" can go
        bounds = self.lp.solve(('C', self.num_cliques),
            bounds=(0, self.num_functions-1))
        return bounds

    def get_clique_bound(self):
        """Gets the bound for finding all the cliques."""
        # ??? omit this?
        bounds = self.get_all_set_bounds()
        return bounds[('C', self.num_cliques)]

if __name__ == '__main__':
    rb = RankBound3(4, 3)
    # lrb.add_average_rank_constraint()
    # print(str(rb.get_all_set_bounds()))
    # print(str(lrb.get_clique_bound()))

