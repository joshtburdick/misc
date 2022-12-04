#!/usr/bin/env python3
# Attempt at bounding the average rank of functions,
# relative to the number of hyperedges.
# ?'s
# - A and C both go "somewhat" high (although not at the top level).
#   On the other hand, B (sets with one edge zonked) are stuck at 0.
#   This seems inaccurate. (For instance, at the very least, we
#   should be able to say A_1==B_1, because they're both the same set!)

import pdb

import scipy.special
import scipy.stats

import lp_helper

class RankBound3:
    """Attempt at bounding the average rank of functions,
    relative to the number of hyperedges.

    """

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        """
        # the number of possible cliques
        self.num_cliques = scipy.special.comb(n, k, exact=True)
        # the number of functions
        self.num_functions = 2 ** self.num_cliques
        # number of cliques which could be "hit" by feeding in a 0 to an edge
        self.max_cliques_zeroed = scipy.special.comb(n-2, k-2, exact=True)
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
            b = scipy.stats.hypergeom(self.num_cliques, i, max_z).pmf(0)
            self.lp.add_constraint(
                [(('C', i), 1.), (('A', i), -(1-b)), (('B', i), -b)],
                '=', 0.)
            # ??? can this be folded into add_zeroing_constraints() ?

    def add_zeroing_constraints(self):
        """Adds constraints from zeroing out one edge."""
        # loop through the "levels" (number of cliques in C)
        # ??? should this start at 1?
        for i in range(self.num_cliques+1):
            # this is the maximum number of cliques which could
            # potentially be zeroed (at this level)
            max_z = min(i, self.max_cliques_zeroed)
            # the probability of some number of cliques being hit
            z = scipy.stats.hypergeom(self.num_cliques, i, max_z)
            # the probability of at least one clique being hit
            p_at_least_one_hit = 1. - z.pmf(0)
            self.lp.add_constraint(
                # The constraint on A (which is "hit") is a weighted sum of
                # what's left in B (none of which are "hit"), after zonking.
                # Note that since this accounting is for A, we normalize
                # by the probability of at least one clique being hit.
                [(('A', i), 1.)] + [
                    (('B', i-j), -z.pmf(j) / p_at_least_one_hit)
                    for j in range(1, max_z+1)],
                '>',
                # Everything in A_i has a higher expected value.
                # Since we're adding distinct things, the expected value
                # of all of the things in A_i is higher
                # (by at least half the number of functions).
                0.5 * p_at_least_one_hit * scipy.special.comb(self.num_cliques, i, exact=True))

    def add_average_rank_constraint(self):
        """Adds equality constraint on average rank of all of the functions."""
        self.lp.add_constraint(
            [(('C', i), scipy.special.comb(self.num_cliques, i, exact=True) / self.num_functions)
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
        bounds = self.lp.solve(('C', self.num_cliques),
            bounds=(0, self.num_functions-1))
        return bounds

    def get_clique_bound(self):
        """Gets the bound for finding all the cliques."""
        # ??? omit this?
        bounds = self.get_all_bounds()
        return bounds[('C', self.num_cliques)]

if __name__ == '__main__':
    rb = RankBound3(7, 4)
    rb.add_mixture_equalities()
    rb.add_zeroing_constraints()
    rb.add_average_rank_constraint()
    print(rb.get_all_bounds())
    # print(str(rb.get_clique_bound()))

