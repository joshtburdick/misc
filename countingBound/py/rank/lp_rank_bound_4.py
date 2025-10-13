#!/usr/bin/env python3
# Attempt at bounding rank, based on zeroing out vertices.
# FIXME
# - use exact rational solver
# - include averages of higher levels
# - use argparse?
# - add sharper upper bounds, based on number of vertices zeroed out?
#   (??? is this still needed? add_vertex_zeroing_constraints()
#   implements this, but doesn't seem to have any effect)

import argparse
import pdb
import sys

import itertools
import numpy as np
import scipy.special
import scipy.stats

sys.path.append("../")            # XXX
sys.path.append("../fractions/")  # XXX

import hypergraph_counter
import lp_helper
import pulp_helper

# Wrapper for comb(), with exact arithmetic.
def comb(n, k):
    return scipy.special.comb(n, k, exact=True)

# Hypergeometric distribution, with exact arithmetic.
def hyperg_frac(N, K, n, k):
    # based on https://en.wikipedia.org/wiki/Hypergeometric_distribution
    # note that we don't try to optimize this
    return fractions.Fraction(
        comb(K, k) * comb(N-K, n-k),
        comb(N, n))

class LpVertexZeroing:
    """Attempt at bound by zeroing out vertices.
    """

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        This will have a variable labelled (v, c) where:
            v is the number of vertices, with k <= v <= n
            c is the number of cliques, with 0 <= c <= {n choose v}
        Each variable (currently) will be the expected rank of functions
        with _up to_ v vertices.

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        # the number of possible cliques
        self.num_cliques = comb(n, k)
        # the total number of functions
        self.num_functions = 2 ** self.num_cliques
        # wrapper for LP solver
        vars = []
        for v in range(k, n+1):
            for c in range(0, comb(v, k) + 1):
                vars += [(v, c)]
        # self.lp = lp_helper.LP_Helper(vars)
        self.lp = pulp_helper.PuLP_Helper(vars)
        # precompute numbers of hypergraphs with different numbers of vertices
        counter = hypergraph_counter.HypergraphCounter(n, k)
        # this is the number which use _up to_ some number of vertices
        self.counts_max_vertices = counter.count_hypergraphs_max_vertices()
        # this is the number which use _all_ of some number of vertices
        self.counts_exact_vertices = counter.count_hypergraphs_exact_vertices()

    def add_average_rank_constraints(self):
        """Adds equality constraints on average rank of all of the functions.

        Idea: possibly this should just be a weighted average of all
        possible relevant functions?
        There is one of these for each possible number of vertices included.

        ??? We could also include something like add_counting_lower_bound(),
        from lp_brute_force_1.py. This changed some of the bounds, but not
        the bound for CLIQUE.
        """
        # add constraints
        for v in range(self.k, self.n+1):
            # get counts of functions, for _up to_ this many vertices
            f = self.counts_max_vertices[v]
            # the number of functions, which implies the bound:
            # both a lower bound (from counting), and an upper bound
            # (because there are vertex-zeroed functions "above" all of these)
            num_functions = f.sum()
            # bound the weighted average of these
            w = f / num_functions
            self.lp.add_constraint(
                [((v, c), w[c]) for c in range(f.shape[0])],
                '=',
                (num_functions-1) / 2.)

    def add_zeroing_upper_bound(self):
        """Adds upper bound based on the number of functions "above".

        ??? This may not be useful; it seems like add_average_rank_constraints()
        almost subsumes it.
        However, this is at least a "hard" bound.
        """
        # loop through number of vertices
        for v in range(self.k, self.n+1):
            # get counts of functions, for _up to_ this many vertices
            f = self.counts_max_vertices[v]
            # the number of functions with up to this many vertices
            num_functions = f.sum()
            # presumably, all of the _other_ functions are "higher"?
            # which (I think) should make this an upper bound;
            # however, this isn't quite clear
            for c in range(0, f.shape[0]):
                # pdb.set_trace()
                self.lp.add_constraint(
                    [((v, c), 1.)],
                    '<=',
                    num_functions)

    def add_vertex_zeroing_constraints(self):
        """Adds constraints from zeroing out vertices.

        The sets of cliques are:
        C: the cliques in the larger set, using (up to) k+1 vertices
        B: the cliques zeroed by feeding in zeros to a vertex of C
        A: the cliques which are left over
        """
        # add case when v == self.k
        # FIXME: this doesn't seem sufficient...
        self.lp.add_constraint([((self.k, 1), 1.)],
            '>=',
            comb(self.n, self.k)/2 + 0.5)

        # loop through number of vertices in larger graph
        for v in range(self.k+1, self.n+1):
            # loop through number of cliques in that graph
            # (note that the trivial bound when C_size==0 is implied by the
            # overall "box" constraints on all the variables)
            for C_size in range(1, comb(v, k)+1):
                # Maximum number of cliques we might hit,
                # _assuming that only v vertices are present_
                # (If we pick a vertex randomly, and miss all of the
                # cliques, we get to "re-roll".)
                num_cliques_hitting_vertex = comb(v-1, k-1)
                # Bounds on number of cliques zeroed.
                # The min is at least 1 (assuming we can "re-roll"), and no
                # more than the difference between the number of cliques in the
                # larger set, and the number left over.
                min_zeroed = max(1, C_size - comb(v-1, k))
                # The max is limited by the current number of vertices (and how
                # many cliques hit one vertex), and the number of cliques.
                max_zeroed = min(num_cliques_hitting_vertex, C_size)
                # the range of possible number of cliques zeroed ...
                B_size = np.arange(min_zeroed, max_zeroed+1)
                # ... and the number left over
                A_size = C_size - B_size
                # the probability of some number of cliques being hit
                # (again, assuming that only v vertices are "in use" by
                # the hyperedges)
                p_zonk = scipy.stats.hypergeom(comb(v, k),
                    C_size,
                    num_cliques_hitting_vertex)
                # the probability of at least one clique being hit
                # (if this happens, we get to "re-roll", so we assume it doesn't)
                # Note that since this we assume at least one clique is hit,
                # we normalize by this probability.
                p_at_least_one_hit = 1. - p_zonk.pmf(0)
                # ??? does this give the same answer?
                # p_at_least_one_hit = np.array([p_zonk.pmf(z1) for z1 in z]).sum()

                # The constraint on the rank of the v-vertex set in C (which is "hit")
                # is a weighted sum of:
                # - the expected rank of what's left in A, after zonking,
                # - half the expected number of functions in A
                # - plus half the number of functions in B
                A = [((v, C_size), 1.)]
                A += [((v-1, C_size-j), -p_zonk.pmf(j) / p_at_least_one_hit)
                    for j in B_size]
                # note that this is the _expected_ number of functions in A
                A_num_functions = ( p_zonk.pmf(B_size) * self.counts_max_vertices[v-1][A_size] ).sum() / p_at_least_one_hit
                # the number of functions (or, "sets of cliques") in B
                # ??? is this right?
                B_num_functions = self.counts_max_vertices[v][ C_size ]
                # pdb.set_trace()
                self.lp.add_constraint(A, '>=',
                    # Since we presumably "hit" a clique, the number of "new"
                    # functions is the number which include all the vertices.
                    (A_num_functions + B_num_functions) / 2.)
 

    def get_bounds_1(self):
        """Gets bounds, with all of the vertices.

        Deprecated (this seems to be very under-determined, and
        it's not clear what these numbers mean).

        Note that this is only minimizing the bound for any number
        of cliques, and all of the vertices.
        However, it returns bounds for any number of cliques or vertices
        (even though it's not explicitly minimizing those individually).
        """
        s = (self.n, self.num_cliques)
        x = self.lp.solve(s, bounds=(0, self.num_functions-1))
        # get bounds for all the vertices, and any number of cliques
        b = np.array([x[(self.n, c)]
            for c in range(self.num_cliques+1)])
        return b

    def get_bounds(self):
        """Gets bounds, for each number of vertices.

        This is the bound at each 'level'.

        ??? is there a more efficient way to compute this?
        """
        def get_bound_at_level(i):
            # this is the number we're minimizing
            x = (self.n, i)
            r = self.lp.solve(x, bounds=(0, self.num_functions-1))
            return r[x]
        # get bounds for all the vertices, and any number of cliques
        b = np.array([get_bound_at_level(c)
            for c in range(self.num_cliques+1)])
        return b

# FIXME add minimizing / maximizing individual variables ?

if __name__ == '__main__':
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    bound = LpVertexZeroing(n, k)

    # FIXME add constraints based on options?
    # ??? do we need the "base case" of only one clique?
    bound.add_average_rank_constraints()
    bound.add_zeroing_upper_bound()
    bound.add_vertex_zeroing_constraints()

    # get bound
    b = bound.get_bounds()

    # print all the bounds, for "all the vertices"
    print('Num. cliques    Expected rank')
    for i in range(bound.num_cliques+1):
        print('\t\t'.join([str(i), str('%.2f' % b[i])]))

