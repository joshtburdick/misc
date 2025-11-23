#!/usr/bin/env python3
# LP bound for clique parity.

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
import flexible_lp_helper
import pulp_helper
import scip_helper
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

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        This will have two groups of variables:
        - tuples of the form ("A", c) where:
            - "c" is the number of cliques, with 0 <= c <= {n choose k}
            - Each variable will be the expected number of gates in
            the sets with exactly c cliques.
        - similarly, tuples of the form ("B", c) where:
            - "c" is the number of cliques after a "bounce" down
        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        self.max_gates = max_gates

        # number of possible cliques
        self.max_cliques = comb(n, k)
        # number of cliques which could be "hit" by zeroing an edge
        self.max_cliques_hit = comb(n-2, k-2)
 
        # the variables
        self.expected_num_gates_vars = (
            [("A", c) for c in range(0, self.max_cliques + 1)]
            + [("B", c) for c in range(0, self.max_cliques + 1 - self.max_cliques_hit)]
        )

        # wrapper for LP solver
        self.lp = pulp_helper.PuLP_Helper(
             self.expected_num_gates_vars + self.num_gates_dist_vars)

        self.basis = gate_basis.TwoInputNandBasis()

        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None

    def add_counting_bound(self, include_B_bound=True):
        """Adds counting bounds, based on number of functions.

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

    def add_step_bound(self):
        """Adds 'step' bounds from the random walk.

        """
        # the bounce "down":



        # the bounce "up":


    def add_level_bound(self):
        """Adds constraints, based on combining "levels" of circuits.

        """
        N = self.num_possible_cliques
        # loop through pairs of sizes we could combine
        for i in range(self.num_possible_cliques+1):
            for j in range(i, self.num_possible_cliques+1):
                # loop through the sizes obtainable by XORing them together
                for xor in range(j-i, min(j+i, self.num_possible_cliques)+1, 2):
                    self.lp.add_constraint(
                        [(("A", xor), 1), (("A", i), -1), (("A", j), -1)],
                        "<=",
                        4)

    def add_no_cliques_constraint(self):
        """Adds trivial constraint, on finding no cliques."""
        # FIXME move this to the "levels" bound
        # ... that is, finding zero cliques requires one NAND gate
        self.lp.add_constraint([(("A", 0), 1)], "=", 1)

    def get_all_bounds(self):
        """Gets bounds for each possible number of cliques.

        This is the bounds for each possible number of cliques,
        in the scenario that the number of gates for
        CLIQUE is minimized.
        """
        # solve, minimizing number of gates for CLIQUE
        r = self.lp.solve(("A", self.num_possible_cliques))
        if not r:
            return None
        # for now, we only get bounds for "expected number of gates"
        # for each number of cliques
        n_cliques = range(self.num_possible_cliques+1)
        bounds = [r[("A", num_cliques)]
            for num_cliques in range(self.num_possible_cliques+1)]
        return pandas.DataFrame({
                'Num. vertices': self.n,
                'Num. cliques': n_cliques,
                'Min. gates': bounds})

def get_bounds(n, k, max_gates, constraints_label,
        use_counting_bound, use_combining_bound, use_no_cliques_bound):
    """Gets bounds with some set of constraints.

    n, k, max_gates: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_step_bound, use_level_bound:
        whether to use each of these groups of constraints
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, max_gates={max_gates}, label={constraints_label}]\n')
    bound = LpParity(n, k)
    if use_counting_bound:
        bound.add_counting_bound()
    if use_step_bound:
        bound.add_level_bound()
    if use_level_bound:
        bound.add_level_bound()
    # pdb.set_trace()
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
    parser.add_argument("--write-transition-matrices",
        help="Dump transition matrices to a file, for debugging",
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
        get_bounds(n, k, max_gates, 'Combining', False, True, False),
        get_bounds(n, k, max_gates, 'Counting and combining', True, True, False),
        get_bounds(n, k, max_gates, 'Counting, combining, and no cliques',
            True, True, True)
    ])
    if args.result_file:
        with open(args.result_file, "wt") as f:
            bounds.to_csv(f, index=False)
    else:
        bounds.to_csv(sys.stdout, index=False)

