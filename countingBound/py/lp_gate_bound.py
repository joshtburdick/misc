#!/usr/bin/env python3
# Linear programming bound, counting number of gates.

print('NOTE: this isn\'t correct.')

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


class TwoInputNandBound:
    """
    Bound on E[ number of two-input NAND gates needed ].

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

    This version tries to do bookkeeping by zeroing out a distinguished edge, e.
    ??? should I rename the edge which is zeroed out?
    ??? should this include a bound everywhere that E[# gates] >= 0 ?
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
        # number of cliques which could be zeroed out when edge e is zeroed out
        self.max_cliques_zeroed = comb(n-2, k-2, exact=True)
        # how many cliques could be left over
        self.max_cliques_remaining = self.max_cliques - self.max_cliques_zeroed
        # mapping from tuples (numVertices, numCliques) to
        # variable index in the LP
        self.var_index = {}
        # set up the mapping of variable indices
        # first, indexed by number of cliques (zeroed, remaining)
        for i in range(self.max_cliques_zeroed+1):
            for j in range(self.max_cliques_remaining+1):
                self.var_index[(i,j)] = len(self.var_index)
        # then, indexed by the total number of cliques
        for i in range(self.max_cliques+1):
            self.var_index[('total_cliques',i)] = len(self.var_index)
        # These store the constraints:
        # A: a list of lists of (A,i,j) entries (which go into a sparse matrix)
        # b: a list of numbers
        # the inequalities (note that the LP solver expects upper bounds)
        self.A_ub = []
        self.b_ub = []
        # the equalities, stored similarly
        self.A_eq = []
        self.b_eq = []
        # counting bound (for this number of inputs)
        num_inputs = comb(n, 2, exact=True)
        self.counting_bound = TwoInputNandBound(num_inputs, 10000)

    def add_constraint(self, A, op, b):
        """Adds one row to the constraints.

        A: a list of (index, coefficient) pairs, where "index" is
            a key (of any hashable Python type) in var_index
        op: either '<', '=', or '>': this is the type of constraint
            ??? can we treat '>' the same as '>='?
        b: the corresponding bound
        Side effects: adds the constraint
        """
        # print(str(A) + ' ' + op + ' ' + str(b))
        # converts from "list of coefficient" to a row of A
        def get_coefs(i, negate): 
            # this is arguably pushing limits for a list comprehension...
            return [(-a if negate else a, i, self.var_index[k])
                for (k,a) in A]
        # add whichever kind of constraint
        if op == '<':
            i = len(self.b_ub)
            self.A_ub += get_coefs(i, False)
            self.b_ub += [b]
            return
        if op == '=':
            i = len(self.b_eq)
            self.A_eq += get_coefs(i, False)
            self.b_eq += [b]
            return
        if op == '>':
            i = len(self.b_ub)
            self.A_ub += get_coefs(i, True)
            self.b_ub += [-b]
            return

    def add_total_cliques_equality_constraints(self):
        """Adds constraints for a given total number of cliques.

        For 0 <= m <= N, these define a variable '(total_cliques, m)',
        which is E[ number of gates need to find m cliques ],
        or "the expected number of gates needed at 'level m'".
        It's constrained to equal the weighted average of FIXME describe this.
        """
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

    def add_total_cliques_counting_bound_constraints(self):
        """Adds counting bound, based on total number of cliques.

        For each "level" of "total number of cliques found", this
        simply adds a lower bound, based on the counting bound.

        ??? Another thing to try is just using the counting bound
        on a weighted average of _all_ the levels. This seems
        simpler in some ways.
        """
        for num_cliques in range(self.max_cliques+1):
            A = [(('total_cliques', num_cliques), 1)]
            b = self.counting_bound.expected_gates(
                math.log2(comb(self.max_cliques, num_cliques, exact=True)))
            self.add_constraint(A, '>', b)

    def add_edge_zeroing_constraints(self):
        """Adds constraints based on zeroing out an edge.

        This is just "if there are A cliques not touching e, and B
        cliques touching it, and C = A \cup B, then E[C] > E[A] + 1" ?
        That seems much simpler, and stronger.
        """
        # loop through number of cliques zeroed (if >= 1)
        for i in range(1, self.max_cliques_zeroed+1):
            # loop through the number of cliques "left over"
            for j in range(self.max_cliques_remaining+1):
                A = [((i,j), 1), ((0,j), -1)]
                # since we're measuring by "# gates", this is
                # "this requires at least one more gate, on average".
                self.add_constraint(A, '>', 1)

    def add_upper_bound_constraints(self):
        """Adds an upper bound.

        """
        # loop through # cliques, with 0 <= a < c <= max_cliques
        for a in range(0, self.max_cliques):
            for c in range(a+1, self.max_cliques+1):
                b = c - a
                # Note that this is an _upper_ bound:
                # |\scriptC(C)| <= |\scriptC(A)| + |\scriptC(B)| + 3
                A = [(('total_cliques', a), -1),
                    (('total_cliques', b), -1),
                    (('total_cliques', c), 1)]
                self.add_constraint(A, '<', 3)

    def solve(self):
        """Solves the linear system.
        
        Note that by default, the solver constrains all x >= 0,
        so we don't add that constraint.
        FIXME: add option to minimize finding only some number
            of cliques (rather than all of them) ?
        Returns: the LP result
        """
        # utility to convert entries to a sparse array
        def sparse_array_from_entries(A):
            # gets i'th element of all the tuples
            def ith(i):
                return [a[i] for a in A]
            return scipy.sparse.coo_array( (ith(0), (ith(1), ith(2))) )
        # convert A and b to np.array objects
        A_ub = sparse_array_from_entries(self.A_ub)
        b_ub = np.array(self.b_ub)
        A_eq = sparse_array_from_entries(self.A_eq)
        b_eq = np.array(self.b_eq)

        # the objective function: how low can the rank of finding
        # all the cliques (with that many vertices) be?
        c = np.zeros(len(self.var_index))
        # c[ self.var_index[(self.max_cliques_zeroed, self.max_cliques_remaining)] ] = 1
        c[ self.var_index[('total_cliques', self.max_cliques)] ] = 1
        # solve
        # ??? Is there a way to tell the solver that this is sparse?
        # (It's detecting this, but that throws a warning.)
        r = scipy.optimize.linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq)
        # FIXME deal with this failing
        
        # pdb.set_trace()
        # Reshape into a rectangle. This is admittedly inefficient when we
        # just want the bound for finding all the cliques; but it seems
        # simplest to just return all of this
        x = np.empty( (self.max_cliques_zeroed+1, self.max_cliques_remaining+1) )
        for i in range(self.max_cliques_zeroed+1):
            for j in range(self.max_cliques_remaining+1):
                x[i,j] = r.x[ self.var_index[(i,j)] ]
        # pdb.set_trace()
        return x

    def get_bound(self, include_upper_bound):
        """Gets the bound.

        include_upper_bound: if True, include the upper bound
        Returns: a 2-D NumPy array, with axes "# cliques zeroed"
            and "# cliques remaining"
        """
        self.add_total_cliques_equality_constraints()
        self.add_total_cliques_counting_bound_constraints()
        self.add_edge_zeroing_constraints()
        # possibly include the upper bound
        if include_upper_bound:
            self.add_upper_bound_constraints()
        x = self.solve()
        return x

def plot_bound_surfaces(n, k, output_prefix):
    """Plots the bound surfaces.

    n, k: size of problem
    output_prefix: prefix for output
    Side effects: plots surfaces, in a file prefixed with
        output_prefix
    """
    print('getting bound without upper bound')
    lp = LpBound(n, k)
    Z0 = lp.get_bound(False)
    print('getting bound with upper bound (this will take longer...)')
    lp = LpBound(n, k)
    Z1 = lp.get_bound(True)
    pdb.set_trace()
    print('plotting surfaces')

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
    parser.add_argument('--plot-surface', action='store_true',
        help='plot the bounds as a surface, with and without the upper bound')
    args = parser.parse_args()

    # possibly plot the surfaces
    if args.plot_surface:
        # plot the surfaces
        plot_bound_surfaces(args.n, args.k, 'bound_surfaces_')
        sys.exit(0)

    # otherwise, print the bound
    lp = LpBound(args.n, args.k)
    x = lp.get_bound(args.include_upper_bound)
    bound = x[ lp.max_cliques_zeroed, lp.max_cliques_remaining ]
    print(np.round(bound, 4))

