#!/usr/bin/env python3
# Attempt based on zeroing out edges (or vertices).

import pdb

import scipy.special
import scipy.stats

import lp_helper


class LpBruteForce1:
    """Attempt at bounding the average rank of functions,

    """

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        This will have one variable for each possible number of cliques
        (or "level"). By convention, level k-1 is the empty set
        of cliques (since with only k-1 vertices, there are no
        k-cliques).

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
        self.lp = lp_helper.LP_Helper(range(k-1, self.num_cliques+1))
        # numbering for the cliques
        cliques = itertools.combinations(range(n), k)
        self.clique_index = dict(zip(cliques, range(len(cliques))))
        # the sets of cliques
        S = np.arange(self.num_functions)
        # the sets of cliques which are "zeroable", and don't overlap S
        Z = np.zeros(self.num_functions)








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
    rb = RankBound3(4, 3)
    # rb.add_zeroing_constraints()
    rb.add_average_rank_constraint()
    pdb.set_trace()
    b = rb.get_all_bounds()
    for i in range(rb.num_cliques+1):
        print('\t'.join([str(i), b[i]))
    print()
    print('bound = ' + str(b[self.num_cliques]))

