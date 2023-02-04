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

import lp_helper

class TwoInputNandBound:
    """
    Bound on E[ num. of 2-input NAND gates needed ], for some num. of functions.

    FIXME factor this out into a separate file?
    """
    def __init__(self, num_inputs, max_log2_num_functions):
        """Constructor (which precomputes a table for this).

        num_inputs: the number of inputs
        Returns a list of tuples of the form (a, b, num_gates), where:
            a, b: are [a, b) interval of a log2 number of functions
            num_gates: number of gates which suffice to any (log_2)
            number of functions in that interval.
        """
        self.num_inputs = num_inputs
        # this tracks the range of number of functions expanded, each
        # time we add a gate; it's mostly for debugging
        self.num_functions_per_gate = []
        # num_gates_needed[i] is (a lower bound on) the number of gates
        # needed to implement 2^i different functions. (It's a lower
        # bound because some of the circuits implement the same function,
        # but using more gates.)
        num_gates_needed = np.full([max_log2_num_functions], np.nan)
        num_gates_needed[0] = 0
        # we will fill this in, adding gates
        # we start with no gates
        num_gates = 0
        # ... which can be used to implement 1 function
        log2_num_functions = 0
        # fill in the table, using more and more gates
        while log2_num_functions < max_log2_num_functions:
            # this is the total number of possible inputs to gates
            m = num_inputs + num_gates
            # the new number of functions expressible, is the previous, plus:
            log2_num_functions_1 = (log2_num_functions
                # this is like "choosing any distinct two of an input, a
                # gate output, or a constant 1 (which converts a two-input
                # NAND gate to simply a NOT gate)"
                # ... plus 1, to allow "less than or equal to this number of gates".
                # This is somewhat over-counting (as one can "waste" two NOT gates
                # to implement a function requiring exactly two fewer gates).
                # However, many distinct circuits will end up computing the same
                # function anyway, which seems likely to cause much more
                # overcounting anyway.
                + math.log2(comb(num_inputs + num_gates + 1, 2, exact=True) + 1))
            # which part of the array to fill in: taking the ceiling
            # seems like a safer assumption
            a = math.ceil(log2_num_functions)
            b = min(math.ceil(log2_num_functions_1), max_log2_num_functions)
            self.num_functions_per_gate.append(
                (num_gates, log2_num_functions, log2_num_functions_1))
            # print(str(a) + ' ' + str(b))
            num_gates_needed[a:b] = num_gates
            num_gates += 1
            log2_num_functions = log2_num_functions_1
        # Now that we have the number of gates needed, we get a lower
        # bound on "expected # of gates needed", by weighting it by
        # the number of functions. (Since adding a wire doubles the
        # number of functions, it's possible that just using
        # "num_gates_needed-1" would suffice.)
        self.expected_gates_needed = np.full([max_log2_num_functions], np.nan)
        # we assume there's an "empty circuit", which computes a 0 or 1
        self.expected_gates_needed[0] = 0
        for i in range(1, max_log2_num_functions):
            # this is an exponential tower, so it's _slightly_ more likely
            # that a random function comes from the top level, than any of
            # the lower levels.
            # Here, we make the conservative assumption that it's
            # equally likely that a random function came from the
            # top layer, or one of the lower layers.
            self.expected_gates_needed[i] = (num_gates_needed[i]
                + self.expected_gates_needed[i-1]) / 2.0

    def expected_gates(self, log2_num_functions):
        """
        Lower bound on expected number of gates.

        log2_num_functions: log_2(number of functions)
        Returns: expected number of gates needed
        """
        # this basically just looks in the table
        return self.expected_gates_needed[ math.floor(log2_num_functions) ]

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
        # number of cliques possible
        self.max_cliques = comb(n, k, exact=True)
        # set up the mapping of variable indices
        # first, variables for circuits with a given total number of cliques

        # FIXME ??? I think that 'total_cliques' goes away,
        # because (I think) it makes sense to just track this
        # for each number of vertices

        var_list = [('total_cliques',i)
            for i in range(self.max_cliques+1)]
        # then, add variables for circuits which find cliques which
        # include a subset of vertices:
        # (i, j) is the average size of circuits which find
        # j cliques, among i vertices

        # FIXME ??? I think the "k-1" case isn't needed?

        # (by convention, "k-1" is no cliques)
        for i in range(k-1, n+1):
            for j in range(comb(self.num_cliques), i)):
                var_list.add((i,j))
        # interface to LP solver (with those variables)
        self.lp = lp_helper(var_list)
        # counting bound (for this number of inputs)
        num_inputs = comb(n, 2, exact=True)
        self.counting_bound = TwoInputNandBound(num_inputs, 10000)

    def add_total_cliques_equality_constraints(self):
        """Adds constraints for a given total number of cliques.

        For 0 <= m <= N, these define a variable '(total_cliques, m)',
        which is E[ number of gates need to find m cliques ],
        or "the expected number of gates needed at 'level m'".
        It's constrained to equal the weighted average of finding
        cliques among subsets of the vertices.
        """

        # FIXME I don't think this is needed

        # loop through the number of cliques
        for num_cliques in range(self.max_cliques+1):
            # bounds on number of cliques containing edge e
            # (these won't actually be zeroed)
            min_cliques_zeroed = max(0, num_cliques - self.max_cliques_remaining)
            max_cliques_zeroed = min(num_cliques, self.max_cliques_zeroed)
            # the probability of some number of cliques containing edge e
            h = hypergeom(
                # number of possible cliques
                self.max_cliques,
                # number of those present
                num_cliques,
                # number of cliques which could intersect edge e
                max_cliques_zeroed)
            # here, z is the number of cliques which _do_ intersect edge e
            A = [((z, num_cliques-z), h.pmf(z))
                for z in range(min_cliques_zeroed, max_cliques_zeroed+1)]
            # this is constraining the total number of gates at this "level"
            # to equal the average, weighted by the probability of some
            # number of cliques being zeroed out
            self.add_constraint(A + [(('total_cliques', num_cliques), -1.0)], '=', 0)

    def add_vertex_zeroing_constraints(self):
        """Adds constraints based on zeroing out a vertex.

        This says that if a function finds some set of cliques C on
        i+1 vertices, and we zero out a vertex, we zero out some
        set B of cliques, and are left with a function which finds
        a set A of those cliques (where A is a proper subset of C),
        and has strictly fewer gates.
        """
        # add "base case" constraints on (k-1,0) ?
        self.lp.add_constraint([(k-1, 0)], '>', 0.)   FIXME
        # loop through the number of vertices in C
        for i in range(self.k+1, self.n+1):
            # this is how many cliques _could_ be in C, and A
            max_C_size = comb(i, k, exact=True)
            max_A_size = comb(i-1, k, exact=True)
            # loop through the number of cliques in C
            for C_size in range(max_C_size+1):
                # distribution of how many cliques in C
                # are in A (as opposed to in B)
                h = hypergeom(max_C_size, max_A_size, C_size)
                # the mixture of levels of A which C is larger than
                A_mix = [((i-1, A_size), -h.pmf(A_size))
                    for A_size in range(max_A_size+1)]
                # C has at least one more gate than what's
                # left in A after zeroing out a vertex
                self.lp.add_constraint(
                    [((i, C_size), 1.)] + A_mix, '>', 1.)
                    

FIXME


                # number of possible cliques
                self.max_cliques,
                # number of those present
                num_cliques,
                # number of cliques which could intersect edge e
                max_cliques_zeroed)
            # here, z is the number of cliques which _do_ intersect edge e
            A = [((z, num_cliques-z), h.pmf(z))
                for z in range(min_cliques_zeroed, max_cliques_zeroed+1)]
 

        # loop through number of cliques zeroed (if >= 1)
        for i in range(1, self.max_cliques_zeroed+1):
            # loop through the number of cliques "left over"
            for j in range(self.max_cliques_remaining+1):
                A = [((i,j), 1), ((0,j), -1)]
                # since we're measuring by "# gates", this is
                # "this requires at least one more gate, on average".
                self.add_constraint(A, '>', 1)

    def add_total_cliques_counting_bound_constraints(self):
        """Adds counting bound, based on total number of cliques.

        For each "level" of "total number of cliques found", this
        simply adds a lower bound, based on the counting bound.

        ??? Another thing to try is just using the counting bound
        on a weighted average of _all_ the levels. This seems
        more useful in some ways.
        """
        for num_cliques in range(self.max_cliques+1):
            A = [(('total_cliques', num_cliques), 1)]
            b = self.counting_bound.expected_gates(
                math.log2(comb(self.max_cliques, num_cliques, exact=True)))
            self.add_constraint(A, '>', b)

    def add_upper_bound_constraints(self):
        """Adds an upper bound.

        This is says that if A and B are sets of cliques, and you have
        circuits to detect each of them, you can detect all of the cliques
        in $A \cup B$, by ORing those circuits together.
        """
        # FIXME add this for each total number of vertices separately
        # loop through # cliques, with 0 <= a < c <= max_cliques
        for a in range(0, self.max_cliques):
            for c in range(a+1, self.max_cliques+1):
                b = c - a
                # Note that this is an _upper_ bound:
                # |\scriptC(C)| <= |\scriptC(A)| + |\scriptC(B)| + 3
                # (the 3 is to OR the circuits for A and B together)
                A = [(('total_cliques', a), -1),
                    (('total_cliques', b), -1),
                    (('total_cliques', c), 1)]
                self.add_constraint(A, '<', 3)

    def get_bound(self):
        """Gets the bound, by solving the linear system.
        
        Note that by default, the solver constrains all x >= 0,
        so we don't add that constraint.
        Returns: the LP result
        """
        # solve
        # ??? Is there a way to tell the solver that this is sparse?
        # (It's detecting this, but throws a warning.)
        x = self.lp.solve(('total_cliques', self.max_cliques))
        # FIXME deal with this failing
        return x[('total_cliques', self.max_cliques)]

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

    lp = LpBound(args.n, args.k)
    self.add_total_cliques_equality_constraints()
    self.add_vertex_zeroing_constraints()
    self.add_total_cliques_counting_bound_constraints()
    # possibly include the upper bound
    if include_upper_bound:
        self.add_upper_bound_constraints()
    # otherwise, print the bound
    bound = lp.get_bound()
    print(np.round(bound, 4))

