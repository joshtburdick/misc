#!/usr/bin/env python3
# Attempt based on zeroing out edges (or vertices).

import argparse
import pdb
import sys

import itertools
import numpy as np
import scipy.special
import scipy.stats

import lp_helper
import hypergraph_counter

def comb(n, k):
    """Short name for comb().

    Note that when n < k, this will return 1."""
    return scipy.special.comb(n, k, exact=True)

class LpVertexZeroing:
    """Attempt at bound by zeroing out vertices.

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
        if n < k:
            raise ValueError('n must be < k')
        # the number of possible cliques
        self.num_cliques = scipy.special.comb(n, k, exact=True)
        # the number of functions
        self.num_functions = 2 ** self.num_cliques
        # wrapper for LP solver, with one variable per "level"
        self.lp = lp_helper.LP_Helper(range(self.num_cliques+1))
        # counts of number of cliques with _exactly_ some number of vertices
        counter = hypergraph_counter.HypergraphCounter(n, k)
        self.hypergraphs_exact_vertices = counter.count_hypergraphs_exact_vertices()

    def add_upper_bound_constraints(self):
        """Adds constraints based on upper bounds from zeroing."""

        # compute upper bound, for each possible number of vertices:
        # any cliques which hit an "unused" vertex are "above" the
        # cliques which only use v vertices
        upper_bound = np.array([
            self.num_functions - (2 ** (self.num_cliques - comb(v, self.k)))
            for v in range(self.n+1)])
        # loop through the "levels" (number of cliques in C)
        for i in range(self.num_cliques+1):
            # number of functions "at this level" (with this many cliques)
            functions_at_level = comb(self.num_cliques, i)
            # get weighted average of upper bound at this level
            average_upper_bound = 0.
            for v in range(self.k-1, self.n+1):
                num_hypergraphs = self.hypergraphs_exact_vertices[v]
                if num_hypergraphs.shape[0] > i:
                    average_upper_bound += (num_hypergraphs[i] / functions_at_level) * upper_bound[v]
            # add constraint at this level
            self.lp.add_constraint([(i, 1.)],
                '<',
                # ... has expected rank less than the upper bound,
                # minus half the number of sets of this size
                average_upper_bound - ((functions_at_level-1) / 2))

    def add_average_rank_constraint(self):
        """Adds equality constraint on average rank of all of the functions."""
        self.lp.add_constraint(
            [(i, scipy.special.comb(self.num_cliques, i, exact=False) / self.num_functions)
                for i in range(self.num_cliques+1)],
            '=',
            (self.num_functions-1) / 2.)

    def add_counting_lower_bounds(self):
        """Adds counting lower bound on a set of lower levels.

        It's not clear whether it's more useful to bound each level i
        individually, or the average of the first i levels.
        (Possibly they work out to the same thing?)
        """
        # this version bounds each level individually
        for i in range(self.num_cliques+1):
            self.lp.add_constraint([(i, 1.)],
                '>',
                (comb(self.num_cliques, i) - 1) / 2)

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

    def get_average_bound_at_top(self, num_top_levels):
        """Gets bound for the top levels.

        Returns: bound for the average of the sets in the
            top num_top_levels levels of sets. (Thus, when
            num_top_levels==1, this is "all the cliques".)
        """
        c = np.zeros(len(self.lp.var_index))
        total_sets = 0.
        for i in range(self.num_cliques+1-num_top_levels, self.num_cliques+1):
            # print(i)
            num_sets_at_level = scipy.special.comb(
                self.num_cliques, i, exact=True)
            c[ i ] = num_sets_at_level
            total_sets += num_sets_at_level
        c = c / total_sets
        x = self.lp.solve_1(c, bounds=(0, self.num_functions-1))
        return {'x':x, 'objective': np.sum(c*x) }

if __name__ == '__main__':
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    # this used to be an arg; for now, just leaving it at 1,
    # to try to minimize rank of CLIQUE
    # num_top_levels = int(sys.argv[3])
    num_top_levels = 1
    bound = LpVertexZeroing(n, k)
    bound.add_upper_bound_constraints()
    bound.add_average_rank_constraint()
    bound.add_counting_lower_bounds()
    b = bound.get_average_bound_at_top(num_top_levels)
    # pdb.set_trace()
    print('Num. cliques    Expected rank')
    for i in range(bound.num_cliques+1):
        print('\t\t'.join([str(i), str('%.2f' % b['x'][i])]))
    #    print(f'{i:5} {.:{20}.{15}}')
    # for now, omitting the average
    # print()
    # print('bound of levels = ' + str(b['objective']))

