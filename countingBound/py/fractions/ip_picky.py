#!/usr/bin/env python3
# Revised IP bound, using the "picky" version of BUGGYCLIQUE.

# FIXME
# - index by num. "yes" and num. "no", rather than num. "total" and num. "no"?
# - add "less efficient" bounds?
# - possibly add_level_constraints() should only include the total number of functions,
#   and add_counting_bound() should be the number of functions which are present?
#   (Mostly, both will be used anyway, so it doesn't necessarily matter that much.)

import argparse
import fractions
import itertools
import math
import pdb
import sys

import numpy as np
import pandas
import scipy.special
import scipy.stats

import gate_basis
import pulp_helper

# Wrapper for comb(), with exact arithmetic.
def comb(n, k):
    return scipy.special.comb(n, k, exact=True)

# Hypergeometric distribution, returning an exact fractions.Fraction .
def hyperg_frac(N, K, n, k):
    # based on https://en.wikipedia.org/wiki/Hypergeometric_distribution
    # note that we don't try to optimize this
    return fractions.Fraction(
        comb(K, k) * comb(N-K, n-k),
        comb(N, n))

class LpPicky:
    """Attempt at bound using PICKYCLIQUE.

    Note that when A and B are sets (of cliques), we write:
    - "A + B" for A union B, and
    - "A <= B" for "A is a (not necessarily proper) subset of B"
    """

    def __init__(self, n, k, max_gates):
        """Constructor gets graph info, and sets up variable names.

        This will have two groups of variables:
        - tuples of the form `("E", c_yes, c_no)` where:
            - c_yes is the number of cliques which must be present for a 1;
              1 <= c_yes <= N (note that we require at least 1 clique to be a "yes")
            - c_no is the number of cliques which, if present, force a 0 output;
              0 <= c_no < N - c_yes
            - Each variable will be the expected number of gates in
              the sets of size `c_yes` and `c_no`.
        - tuples of the form `(c_yes, c_no, g)`, where
            - c_yes and c_no are the numbers of cliques (as above)
            - g is the number of gates
            - Each variable will be the number of circuits having g gates, which
              detect c_yes cliques, AND NOT c_no cliques
        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        max_gates: maximum number of gates to include
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        self.max_gates = max_gates
        # number of possible cliques
        self.num_possible_cliques = comb(n, k)

        self.expected_num_gates_vars = []
        self.num_gates_dist_vars = []
        # number of "yes" cliques (at least one of which must be present
        # for a 1 to be output)
        for i in range(1, self.num_possible_cliques + 1):
            # number of "no" cliques which, if present, force a 0 output
            # (even if a "yes" clique is present)
            for j in range( self.num_possible_cliques + 1 - i ):
                # variables for expected number of gates, for each number of cliques
                self.expected_num_gates_vars += [("E", i, j)]
                for g in range(1, max_gates+1):
                    # variables for counts of numbers of functions
                    # with some number of gates
                    self.num_gates_dist_vars += [(i, j, g)]

        # wrapper for LP solver
        # FIXME: make this an option?
        self.lp = pulp_helper.PuLP_Helper(
             self.expected_num_gates_vars + self.num_gates_dist_vars)

        # FIXME make this an option?
        self.basis = gate_basis.UnboundedFanInNandBasis()
        # self.basis = gate_basis.TwoInputNandBasis()

        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None

    def add_level_constraints(self):
        """Adds constraints on functions at some "level".

        By "level", we mean "number of cliques".

        The constraints (both of which are equalities) are:
        - on the total number of functions at that "level", and
        - connecting the counts with that "level"'s expected gate count
        """
        # loop through number of cliques
        for i in range(self.num_possible_cliques+1):
            for j in range(i-1):
                # we compute this by:
                # - choosing a set of "yes" cliques, then
                # - choosing a set of "no" cliques, from those which remain
                num_functions = (comb(self.num_possible_cliques, i)
                    * comb(self.num_possible_cliques - i, j))
                # add constraint that these sum to the number of functions
                self.lp.add_constraint(
                    [((i, j, g), 1) for g in range(1, self.max_gates+1)],
                    '=', num_functions)
                # add constraint defining expected number of gates
                A = [((i, j, g), g)
                    for g in range(1, self.max_gates+1)]
                self.lp.add_constraint(A + [(("E", i, j), -num_functions)],
                    '=', 0)

    def add_counting_bound(self):
        """Adds counting bounds, for a given number of gates.

        This is a lower bound on the number of functions with
        some number of gates.
        """
        # number of possible functions for each possible number of gates
        # (with number of inputs based on number of vertices)
        num_possible_functions = self.basis.num_functions(comb(self.n, 2), self.max_gates+1)
        # upper-bound "total number of functions with this many gates"
        for g in range(1, self.max_gates+1):
            # this list comprehension is admittedly baroque
            self.lp.add_constraint(
                [((i, j, g), 1)
                    for i in range(1, self.num_possible_cliques + 1)
                    for j in range(self.num_possible_cliques + 1 - i)],
                '<=', num_possible_functions[g])

    def add_buggy_bound(self):
        """Adds upper bound on computing 'buggy' sets of functions.

        We can implement BUGGYCLIQUE(A+B) by combining circuits:
        PICKYCLIQUE(A,B) OR BUGGYCLIQUE(D), where B <= D <= A+B.
        """
        for i in range(1, self.num_possible_cliques + 1):
            for j in range(self.num_possible_cliques + 1 - i):
                for k in range(j, i+j+1):
                    self.lp.add_constraint(
                        [(("E",i+j,0), 1), (("E",i,j), -1), (("E",k,0), -1)],
                        "<=",
                        # Note that for unbounded fan-in, "OR" can be implemented
                        # by combining all of the inputs to the last gate in
                        # each circuit.
                        -1)    # FIXME the "OR" cost should depend on the basis

    def add_picky_bound(self):
        """Adds upper bound on computing 'picky' sets of functions.

        We can implement PICKYCLIQUE(A,B) by combining circuits:
        BUGGYCLIQUE(D) AND NOT BUGGYCLIQUE(B), where A <= D <= A+B.
        """
        for i in range(1, self.num_possible_cliques + 1):
            for j in range(self.num_possible_cliques + 1 - i):
                for k in range(i, i+j+1):
                    self.lp.add_constraint(
                        [(("E",i,j), 1), (("E",k,0), -1), (("E",j,0), -1)],
                        "<=",
                        3)    # FIXME the "AND NOT" cost should depend on the basis
                        # FIXME double check "AND NOT" cost?

    def add_naive_upper_bound(self):
        """Adds naive upper bound."""
        # add bound for finding one or more cliques
        for num_cliques in range(1, self.num_possible_cliques+1):
            # FIXME should depend on basis
            num_gates = num_cliques + 1
            if num_gates <= self.max_gates:
                self.lp.add_constraint([(("E",num_cliques,0), 1)], "<=", num_gates)
            else:
                # stop when we hit the maximum number of gates being considered
                return

    def get_all_bounds(self):
        """Gets bounds for each possible number of cliques.

        This is the bounds for each possible number of cliques,
        in the scenario that the number of gates for
        CLIQUE is minimized.
        """
        # solve, minimizing number of gates for CLIQUE
        r = self.lp.solve(("E", self.num_possible_cliques, 0))
        if not r:
            return None
        # XXX for now, just getting counts for BUGGYCLIQUE
        # (with 0 'no's)
        n_cliques = range(self.num_possible_cliques+1)
        bounds = [r[("E", num_cliques, 0)]
            for num_cliques in range(1, self.num_possible_cliques+1)]
        return pandas.DataFrame({
            'Num. vertices': self.n,
            'Num. cliques': n_cliques,
            'Min. gates': bounds})

def get_bounds(n, k, max_gates, constraints_label,
        use_counting_bound, use_picky_bound, use_upper_bound):
    """Gets bounds with some set of constraints.

    n, k, max_gates: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_buggy_bound, use_picky_bound, use_upper_bound:
        whether to use each of these groups of constraints
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, max_gates={max_gates}, label={constraints_label}]\n')
    bound = LpPicky(n, k, max_gates)
    bound.add_level_constraints()
    if use_counting_bound:
        bound.add_counting_bound()
    if use_buggy_bound:
        bound.add_picky_bound()
    if use_picky_bound:
        bound.add_picky_bound()
    if use_upper_bound:
        bound.add_naive_upper_bound()
    b = bound.get_all_bounds()
    b['Constraints'] = constraints_label
    return b.iloc[:,[3,0,1,2]]

def parse_args():
    parser = argparse.ArgumentParser(
        description='Bounds circuit size for detecting subsets of cliques.')
    parser.add_argument("n", type=int,
        help="Number of vertices")
    parser.add_argument("k", type=int,
        help="Size of clique")
    parser.add_argument("max_gates", type=int,
        help="Maximum number of gates to consider (if this is too"
            " small, thr problem will be infeasible)")
    parser.add_argument("--dump-lp",
        help="Dump LP problem statement to a file",
        action="store_true")
    parser.add_argument("--result-file",
        help="Write result to indicated file (rather than stdout)")
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    n = args.n
    k = args.k
    max_gates = args.max_gates

    # output_file = sys.argv[3]
    bounds = pandas.concat([
        get_bounds(n, k, max_gates, 'Counting', True, False, False, False),
        get_bounds(n, k, max_gates, 'Counting, buggy, and picky', True, True, True, False),
        get_bounds(n, k, max_gates, 'Counting, buggy, picky, and upper bound', True, True, True, True),
    ])
    if args.result_file:
        with open(args.result_file, "wt") as f:
            bounds.to_csv(f, index=False)
    else:
        bounds.to_csv(sys.stdout, index=False)

