#!/usr/bin/env python3
# IP bound for clique parity: simplified.

import argparse
import fractions
import itertools
import math
import pdb
import sys

sys.path.append("..")   # XXX

import numpy as np
import pandas
import scipy.special
import scipy.stats

import gate_basis
import flexible_lp_helper
import pulp_helper
import scip_helper
# import exact_simplex_helper
import simplex_algorithm_helper

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

class LpParity:
    """Attempt at bound for the clique parity problem.
    """

    def __init__(self, n, k, max_gates):
        """Constructor gets graph info, and sets up variable names.

        This will have two groups of variables:
        - tuples of the form ("E", c) where:
            - "c" is the number of cliques, with 0 <= c <= {n choose k}
            - Each variable will be the expected number of gates in
            the sets with exactly c cliques.
        - tuples of the form (c, g), where
            - "c" is the number of cliques
            - "g" is the number of gates
            - Each variable will be the number of sets of "c" cliques,
              having "g" gates.
        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        max_gates: maximum number of gates to include
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        self.max_gates = max_gates

        self.expected_num_gates_vars = []
        self.num_gates_dist_vars = []
        for c in range(0, comb(n, k) + 1):
            # variables for expected number of gates, for each number of cliques
            self.expected_num_gates_vars += [("E", c)]
            for g in range(1, max_gates+1):
                # variables for counts of numbers of functions
                # with some number of gates
                self.num_gates_dist_vars += [(c, g)]
        # wrapper for LP solver
        self.lp = pulp_helper.PuLP_Helper(
             self.expected_num_gates_vars + self.num_gates_dist_vars)
        # self.lp = simplex_algorithm_helper.SimplexAlgorithmHelper(
        #      self.expected_num_gates_vars + self.num_gates_dist_vars, verbosity=0)

        # number of possible cliques
        self.num_possible_cliques = comb(n, k)

        self.basis = gate_basis.TwoInputNandBasis()
        self.rng = np.random.default_rng()

        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None

    def add_level_constraints(self):
        """Adds constraints on functions at some "level"

        By "level", we mean "number of cliques".

        The constraints (both of which are equalities) are:
        - on the total number of functions at that "level", and
        - connecting the counts with that "level"'s expected gate count
        """
        # loop through number of cliques
        for c in range(self.num_possible_cliques+1):
            num_functions = comb(self.num_possible_cliques, c)
            # add constraint that these sum to the number of functions
            self.lp.add_constraint(
                [((c, i), 1) for i in range(1, self.max_gates+1)],
                '=', num_functions)
            # add constraint defining expected number of gates
            A = [((c, i), i)
                for i in range(1, self.max_gates+1)]
            self.lp.add_constraint(A + [(("E", c), -num_functions)],
                '=', 0)

    def add_counting_bound(self):
        """Adds counting bounds, for a given number of gates.

        This is a lower bound on the number of functions with
        some number of gates.
        """
        # number of possible functions for each possible number of gates
        # (with number of inputs based on number of vertices)
        num_functions = self.basis.num_functions(comb(self.n, 2), self.max_gates+1)
        # upper-bound "total number of functions with this many gates"
        for g in range(1, self.max_gates+1):
            self.lp.add_constraint(
                [((c, g), 1)
                    for c in range(self.num_possible_cliques+1)],
                '<=', num_functions[g])

    def num_gates_upper_bound(self, num_cliques):
        """Naive upper bound on computing parity of some number of cliques."""
        if num_cliques == 0:
            return 1
        # number of gates to detect each individual clique, by ANDing the edges
        gates_for_cliques = num_cliques * 2 * (comb(k,2)-1)
        # number of gates to XOR these together
        gates_for_xor = 4 * (num_cliques-1)
        return gates_for_cliques + gates_for_xor

    def add_small_constraint(self):
        """Adds constraint on computing parity of small sets of cliques."""
        for i in range(self.num_possible_cliques // 2):
            self.lp.add_constraint([(("E", i), 1)],
                "<=",
                self.num_gates_upper_bound(i))

    def add_large_constraint(self):
        """Adds constraint on computing parity of large sets of cliques."""
        N = self.num_possible_cliques
        for i in range(1, N // 2):
            # this bounds "the number of gates to find parity of 'a few less'
            # than N cliques", relative to finding the parity of all of them
            self.lp.add_constraint([(("E", N-i), 1), (("E", N), -1)],
                "<=", 
                self.num_gates_upper_bound(i) + 4)

    def get_all_bounds(self):
        """Gets bounds for each possible number of cliques.

        This is the bounds for each possible number of cliques,
        in the scenario that the number of gates for
        CLIQUE is minimized.
        """
        # solve, minimizing number of gates for CLIQUE
        r = self.lp.solve(("E", self.num_possible_cliques))
        if not r:
            return None
        # for now, we only get bounds for "expected number of gates"
        # for each number of cliques
        n_cliques = range(self.num_possible_cliques+1)
        bounds = [r[("E", num_cliques)]
            for num_cliques in range(self.num_possible_cliques+1)]
        return pandas.DataFrame({
                'Num. vertices': self.n,
                'Num. cliques': n_cliques,
                'Min. gates': bounds})

def get_bounds(n, k, max_gates, constraints_label,
        use_counting_bound, use_small_bound, use_large_bound):
    """Gets bounds with some set of constraints.

    n, k, max_gates: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_small_bound, use_large_bound:
        whether to use each of these groups of constraints
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, max_gates={max_gates}, label={constraints_label}]\n')
    bound = LpParity(n, k, max_gates)
    bound.add_level_constraints()
    if use_counting_bound:
        bound.add_counting_bound()
    if use_small_bound:
        bound.add_small_constraint()
    if use_large_bound:
        bound.add_large_constraint()
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
        help="Maximum number of gates to consider")
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

    bounds = pandas.concat([
        get_bounds(n, k, max_gates, 'Counting', True, False, False),
        get_bounds(n, k, max_gates, 'Counting and small', True, True, False),
        get_bounds(n, k, max_gates, 'Counting and large', True, False, True),
        get_bounds(n, k, max_gates, 'Counting, small and large',
            True, True, True)
    ])
    if args.result_file:
        with open(args.result_file, "wt") as f:
            bounds.to_csv(f, index=False)
    else:
        bounds.to_csv(sys.stdout, index=False)

