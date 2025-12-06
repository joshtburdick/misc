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

# Binomial distribution, returning an exact fractions.Fraction .
def binom_frac(n, k):
    return fractions.Fraction(comb(n, k), 2**n)

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

        # number of inputs to the circuit
        self.num_inputs = comb(n, 2)
        # number of possible cliques
        self.max_cliques = comb(n, k)
        # number of cliques which could be "hit" by zeroing an edge
        self.max_cliques_hit = comb(n-2, k-2)
        # number of edges per clique
        self.edges_per_clique = comb(k, 2)
        # the variables
        self.expected_num_gates_vars = (
            [("A", c) for c in range(0, self.max_cliques + 1)]
            + [("B", c) for c in range(0, self.max_cliques + 1 - self.max_cliques_hit)]
        )

        # wrapper for LP solver
        self.lp = pulp_helper.PuLP_Helper(
             self.expected_num_gates_vars)

        self.basis = gate_basis.TwoInputNandBasis()

        # Upper bound on number of gates hit by zeroing an edge
        # (as a function of number of cliques hit).
        self.num_gates_hit_upper_bound = np.array([0]
            + [(self.basis.xor_of_and_upper_bound(self.edges_per_clique, i)
                + self.basis.xor_upper_bound())
                for i in range(1, self.max_cliques_hit + 1)])

        # Lower bound on number of gates hit by zeroing an edge
        # (again, as a function of number of cliques hit).        
        self.num_gates_hit_lower_bound = np.array([0] + [1] * self.max_cliques_hit)

        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None

    def add_counting_bound(self, include_B_bound=True):
        """Adds counting bounds, based on number of functions.

        For each "level" of number of cliques, we bound the
        number of functions with that many cliques.
        """
        # first, compute the number of gates for however many functions have
        # a given number of cliques
        num_functions = np.array([comb(self.max_cliques, i)
                for i in np.arange(self.max_cliques+1)])
        num_gates = self.basis.expected_num_gates(self.num_inputs,
            num_functions)
        for i in range(self.max_cliques+1):
            self.lp.add_constraint(
                [(("A", i), 1)],
                '>=',
                num_gates[i])
        # similar bound, for B variables (after bounce down)
        num_gates = self.basis.expected_num_gates(self.num_inputs,
            np.array([comb(self.max_cliques+1-self.max_cliques_hit, i)
                for i in np.arange(self.max_cliques+1-self.max_cliques_hit)]))
        for i in range(self.max_cliques+1-self.max_cliques_hit):
            self.lp.add_constraint(
                [(("B", i), 1)],
                '>=',
                num_gates[i])

    def add_bounce_down_bound(self):
        """Adds 'bounce down' bound.

        This bounds the number of gates for CLIQUE-PARITY
        for different numbers of cliques, after one edge
        is zeroed out.
        """
        # loop through possible numbers of cliques after bounce down
        for i in range(self.max_cliques-self.max_cliques_hit+1):
            # First, compute the binomial distribution, for the number
            # of cliques before the bounce down.
            p = np.array([binom_frac(self.max_cliques_hit, j)
                for j in range(self.max_cliques_hit+1)])
            # the coefficients for the constraint
            coefs = [(("B", i), 1)] + [(("A", i+j), -p[j])
                for j in range(p.shape[0])]
            # The lower bound: we can't have lost more gates than
            # (roughly) a constant times the number of cliques we hit.
            self.lp.add_constraint(
                    coefs,
                    ">=",
                    (-p * self.num_gates_hit_upper_bound).sum())
            # the upper bound: we hit at least one gate.
            self.lp.add_constraint(
                    coefs,
                    "<=",
                    (-p * self.num_gates_hit_lower_bound).sum())

    def add_bounce_up_bound(self):
        """Adds 'bounce up' bound.

        """
        # loop through possible numbers of cliques after bounce up
        for i in range(self.max_cliques+1):
            # First, compute the range of how many cliques we could
            # have added on the bounce up.
            j_range = np.arange(max(0, i-(self.max_cliques-self.max_cliques_hit)),
                            min(i, self.max_cliques_hit)+1)
            # Next, compute the hypergeometric distribution for how likely each of those is.
            p = np.array([hyperg_frac(self.max_cliques, self.max_cliques_hit, i, j)
                for j in j_range])
            print(f"i = {i}, j_range = {j_range}, p = {p}")
            # Now, add the constraints.
            coefs = [(("A", i), 1)] + [(("B", i-j_range[k]), -p[k])
                for k in range(p.shape[0])]
            # The lower bound: we (probably) added at least one gate.
            self.lp.add_constraint(
                coefs,
                ">=",
                (p * self.num_gates_hit_lower_bound[j_range]).sum())
            # The upper bound: we can't have added more gates than
            # (roughly) a constant times the number of cliques we added.
            self.lp.add_constraint(
                coefs,
                "<=",
                (p * self.num_gates_hit_upper_bound[j_range]).sum())


    def add_combining_bound(self):
        """Adds constraints, based on combining "levels" of circuits.

        """
        # loop through pairs of sizes we could combine
        for i in range(self.max_cliques):
            for j in range(i+1, self.max_cliques+1):
                # loop through the sizes obtainable by XORing them together
                for xor in range(j-i, min(j+i, self.max_cliques)+1, 2):
                    self.lp.add_constraint(
                        [(("A", xor), 1), (("A", i), -1), (("A", j), -1)],
                        "<=",
                        # We need to XOR these together.
                        self.basis.xor_upper_bound)

    def add_no_cliques_constraint(self):
        """Adds trivial constraint, on finding no cliques."""
        # FIXME move this to the "levels" bound
        # ... that is, finding zero cliques requires no NAND gate
        self.lp.add_constraint([(("A", 0), 1)], "=", 0)

    def get_all_bounds(self):
        """Gets bounds for each possible number of cliques.

        This is the bounds for each possible number of cliques,
        in the scenario that the number of gates for
        CLIQUE is minimized.
        """
        # solve, minimizing number of gates for CLIQUE
        r = self.lp.solve(("A", self.max_cliques))
        if not r:
            return None
        # for now, we only get bounds for "expected number of gates"
        # for each number of cliques
        n_cliques = range(self.max_cliques+1)
        bounds = [r[("A", num_cliques)]
            for num_cliques in range(self.max_cliques+1)]
        return pandas.DataFrame({
                'Num. vertices': self.n,
                'Num. cliques': n_cliques,
                'Min. gates': bounds})

def get_bounds(n, k, constraints_label,
        use_counting_bound, use_step_bound, use_combining_bound):
    """Gets bounds with some set of constraints.

    n, k: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_step_bound, use_level_bound:
        whether to use each of these groups of constraints
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, label={constraints_label}]\n')
    bound = LpParity(n, k)

    if use_counting_bound:
        bound.add_counting_bound()
    if use_step_bound:
        bound.add_bounce_down_bound()
        bound.add_bounce_up_bound()
    if use_combining_bound:
        bound.add_combining_bound()
    b = bound.get_all_bounds()
    if b is None:
        return None
    b['Constraints'] = constraints_label
    return b.iloc[:,[3,0,1,2]]

def parse_args():
    parser = argparse.ArgumentParser(
        description='Bounds circuit size for detecting subsets of cliques.')
    parser.add_argument("n", type=int,
        help="Number of vertices")
    parser.add_argument("k", type=int,
        help="Size of clique")
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
    bounds = pandas.concat([
        get_bounds(n, k, 'Counting', True, False, False),
        get_bounds(n, k, 'Counting and step', True, True, False),
        # get_bounds(n, k, 'Counting, step, and combining', True, True, True),
    ])
    if args.result_file:
        with open(args.result_file, "wt") as f:
            bounds.to_csv(f, index=False)
    else:
        bounds.to_csv(sys.stdout, index=False)

