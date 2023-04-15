#!/usr/bin/env python3
# Linear programming bound, counting number of gates.

import argparse
import itertools
import math
import pdb
import sys

import numpy as np
import scipy.optimize
import scipy.sparse

# note that comb() returns a float by default;
# for loop bounds, it needs the "exact=True" option,
# so that it returns an int
from scipy.special import comb
from scipy.stats import hypergeom

import function_gate_count
import lp_helper
import hypergraph_counter

def count_bits(x):
    """Counts number of bits set in a numpy vector."""
    # we assume numbers are "somewhat small" (if they were
    # large, they'd take a while to loop through anyway)
    num_bits_set = np.zeros(x.shape[0])
    for i in range(30):
        mask = np.array(2**i)
        num_bits_set += ((x & mask) > 0) + 0
    return num_bits_set

class LpBound:
    """Computes a bound on the number of gates for finding cliques.

    Here, we zero out a vertex at a time (as zeroing out edges seems
    complicated).
    """
    def __init__(self, n, k):
        """ Constructor.

        ??? should I rename these?
        n: number of vertices
        k: size of clique we're looking for
        """
        # problem size
        self.n = n
        self.k = k
        # set up the mapping of variable indices
        var_list = []
        # add variables for circuits which find cliques which
        # include some number of vertices:
        # (i, j) is the average size of circuits which find
        # j cliques, among i vertices
        for i in range(k, n+1):
            num_cliques = comb(i, self.k, exact=True)
            for j in range(comb(i, k, exact=True)+1):
                var_list.append((i,j))
        # interface to LP solver (with those variables)
        self.lp = lp_helper.LP_Helper(var_list)

    def add_vertex_zeroing_constraints(self):
        """Adds constraints based on zeroing out a vertex.

        This says that if a function finds some set of cliques C on
        i+1 vertices, and we zero out a vertex, we zero out some
        set B of cliques, and are left with a function which finds
        a set A of those cliques (where A is a proper subset of C),
        and has strictly fewer gates.
        """
        # add "base case" constraints on k vertices
        self.lp.add_constraint([((self.k, 0), 1.)], '>', 0.)
        # ??? use the bound that this is at least $(k-1) \choose 2$,
        # just because all the inputs need to be used?
        # it doesn't seem like it would help much.
        self.lp.add_constraint([((self.k, 1), 1.)], '>', 1.)

        # loop through the number of vertices in C
        # (which has at least k+1 vertices)
        for i in range(self.k+1, self.n+1):
            # this is how many cliques _could_ be in C and A, respectively
            max_C_size = comb(i, self.k, exact=True)
            max_A_size = comb(i-1, self.k, exact=True)
            # loop through the number of cliques in C
            for C_size in range(max_C_size+1):
                # probability of a given number of cliques
                # being in A (as opposed to B)
                h = hypergeom(max_C_size, max_A_size, C_size)
                # the mixture of levels of A which C is larger than
                A_mix = [((i-1, A_size), -h.pmf(A_size))
                    for A_size in range(max_A_size+1)]
                # Usually, C has at least one more gate than what's
                # left in A after zeroing out a vertex. However, if
                # nothing is zonked (so |B| == 0, and |A| == |C|),
                # then we're only guaranteed that C requires at least
                # as many gates as A. Thus, for now, we only encode
                # that C's circuit is at least as large.
                self.lp.add_constraint(
                    [((i, C_size), 1.)] + A_mix, '>', 0.)

    def add_total_cliques_counting_bound_constraints(self):
        """Adds counting bound, based on total number of functions.

        For each number of vertices, and for each "level" of
        "total number of cliques found", this simply adds a lower bound,
        based on the counting bound.

        ??? Another thing to try is just using the counting bound on:
        - a weighted average of _all_ the levels, or
        - a weighted average of the levels _up to_ some size
        Maybe this is useful.

        ??? Can we count the fact that we can zero out any of the
        vertices? (In other words, are the vertices "labeled"?)
        It seems like we should be able to: if you're finding
        triangles in a 6-vertex graph, you can zonk distinct
        sets of vertices, and end up with distinct triangles.
        However, the book-keeping seems more complicated.
        """
        # loop through number of vertices
        for num_vertices in range(self.k, self.n+1):
            # number of cliques possible, for this many vertices
            max_cliques = comb(num_vertices, self.k, exact=True)

            # ??? should this be the counting bound for num_vertices?
            # or is it for n? for now, this is for num_vertices
            counting_bound = function_gate_count.TwoInputNandBound(
                # the number of input wires
                comb(num_vertices, 2, exact=True),
                # log_2 of the number of functions shouldn't exceed this
                max_cliques + 5)
            # loop through how many cliques could be present
            for num_cliques in range(max_cliques+1):
                A = [((num_vertices, num_cliques), 1.)]
                b = counting_bound.expected_gates(
                    math.log2(comb(max_cliques, num_cliques, exact=True)))
                self.lp.add_constraint(A, '>', b)

    def add_counting_bound_2(self):
        """Adds counting bound, based on total number of functions.

        For each number of vertices, and for each "level" of
        "total number of cliques found", this simply adds a lower bound,
        based on the counting bound.

        This assumes that we can we count the fact that we can zero out any of the
        vertices; in other words, the vertices are "labeled".
        It seems like we should be able to: if you're finding
        triangles in a 6-vertex graph, you can zonk distinct
        sets of vertices, and end up with distinct triangles.

        The book-keeping seems a bit complicated.
        """
        # in this version of the counting bound, the number
        # of input wires is the same (even when some vertices
        # have been zonked).
        counting_bound = function_gate_count.TwoInputNandBound(
            # the number of input wires
            comb(self.n, 2, exact=True),
            # log_2 of the number of functions shouldn't exceed this
            comb(self.n, self.k, exact=True) + 5)

        # count of functions
        hc = hypergraph_counter.HypergraphCounter(self.n, self.k)
        num_functions = hc.count_hypergraphs_max_vertices()

        # loop through number of vertices not "zonked"
        for num_vertices in range(self.k, self.n+1):
            # number of cliques possible, for this many vertices
            max_cliques = comb(num_vertices, self.k, exact=True)

            # loop through how many cliques could be present
            for num_cliques in range(max_cliques+1):
                A = [((num_vertices, num_cliques), 1.)]
                b = counting_bound.expected_gates(
                    math.log2(num_functions[num_vertices][num_cliques]))
                self.lp.add_constraint(A, '>', b)

    def add_upper_bound_constraints(self):
        """Adds an upper bound.

        This says that if A and B are sets of cliques, and you have
        circuits to detect each of them, you can detect all of the cliques
        in $A \cup B$, by ORing those circuits together.

        Again, there's a copy of this for each number of vertices.
        """
        # loop through number of vertices
        for num_vertices in range(self.k, self.n+1):
            max_cliques = comb(num_vertices, self.k, exact=True)
            # loop through # cliques, with 0 <= a < c <= max_cliques
            for a in range(0, max_cliques):
                for c in range(a+1, max_cliques+1):
                    b = c - a
                    # Note that this is an _upper_ bound:
                    # |\scriptC(C)| <= |\scriptC(A)| + |\scriptC(B)| + 3
                    # (the 3 is to OR the circuits for A and B together)
                    A = [((num_vertices, a), -1.),
                        ((num_vertices, b), -1.),
                        ((num_vertices, c), 1.)]
                    self.lp.add_constraint(A, '<', 3)

    def get_bound(self):
        """Gets the bound, by solving the linear system.
        
        Note that by default, the solver constrains all x >= 0,
        so we don't add that constraint.
        Returns: the LP result
        """
        # index of "all the cliques"
        i = (self.n, comb(self.n, self.k, exact=True))
        # solve
        # ??? Is there a way to tell the solver that this is sparse?
        # (It's detecting this, but throws a warning.)
        x = self.lp.solve(i)
        # FIXME deal with this failing
        return x

def gate_bound_smoke_test():
    """Basic test of 2-input NAND gate counting bound.

    FIXME check these numbers
    """
    counting_bound = TwoInputNandBound(3, 60)
    for b in counting_bound.num_functions_per_gate:
        print(b)

if __name__ == '__main__':
    # gate_bound_smoke_test()

    parser = argparse.ArgumentParser(
        description='Attempt at bound on finding some number of cliques.')
    parser.add_argument('n', type=int,
        help='number of vertices in input graph')
    parser.add_argument('k', type=int,
        help='number of vertices in cliques to find')
    parser.add_argument('--include-upper-bound', action='store_true',
        help='include the upper bound constraint')
    args = parser.parse_args()

    bound = LpBound(args.n, args.k)
    bound.add_vertex_zeroing_constraints()

    bound.add_total_cliques_counting_bound_constraints()
    # revised version of this
    # bound.add_counting_bound_2()

    # possibly include the upper bound
    if args.include_upper_bound:
        bound.add_upper_bound_constraints()
    # otherwise, print the bound
    # FIXME print bound for all numbers of cliques (not just all of them)
    x = bound.get_bound()
    for ((n, k), b) in x.items():
        print(f'{n} {k} {np.round(b, 4)}')
    print()
    # key for "all the cliques"
    i = (bound.n, comb(bound.n, bound.k, exact=True))
    print('overall bound = ' + str(np.round(x[i], 4)))

