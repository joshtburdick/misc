#!/usr/bin/env python3
# LP bound for clique parity, with layers.

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

class IpLayers:
    """Attempt at bound for "layers" of the clique problem.
    """

    def __init__(self, n, k, max_gates, num_layers):
        """Constructor gets graph info, and sets up variable names.

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        max_gates: maximum number of gates to consider (if this is too small,
            the LP will be infeasible)
        num_layers: If this is an integer, the number of layers to use.
            For now, this should be odd.
        """
        self.n = n
        self.k = k
        self.max_gates = max_gates
        if k < 3:
            raise ValueError('k must be >= 3')

        # Number of possible cliques
        self.num_possible_cliques = comb(n, k)

        # The layers, as tuples (a, b), which represent the interval [a, b).
        # Note that we require these to be symmetric around N/2.
        self.layers = self.make_layers(num_layers)

        # The variables for the LP.
        vars = []
        # Counts for each number of vertices,
        # layer (denoted by just the lower endpoint), and number of gates.
        for v in range(self.k, self.n+1):
            max_cliques_for_v = comb(v, self.k)
            for l1, l2 in self.layers:
                if l2 <= max_cliques_for_v:
                    vars += [((v, l1, g) for g in range(self.max_gates+1)]
        # Average for each number of vertices.
        for v in range(self.k, self.n+1):
            vars += [("V", v)]
        # Average for each layer.
        for l1, l2 in self.layers:
            vars += [("L", l1)]

        # wrapper for LP solver
        self.lp = pulp_helper.PuLP_Helper(vars)

        self.basis = gate_basis.TwoInputNandBasis()

        # Counts of functions in each layer, with each possible
        # number of vertices.
        


        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None




    def make_layers(self, num_layers):
        """Makes the layers for the LP.
        
        We divide the possible number of cliques into layers, and add
        constraints on the average number of gates for each layer.
        We require the layers to be symmetric around N/2.
        (The layers won't, in general, be the same height; but this
        tries to make them as similar as possible.)

        num_layers: currently must be odd.
        """
        if num_layers % 2 == 0:
            raise ValueError('num_layers must be odd')
        height = self.max_cliques // num_layers
        endpoints = range(0, self.max_cliques, height)
        endpoints += [self.max_cliques-i for i in reversed(endpoints)]
        return [(endpoints[i], endpoints[i+1]) for i in range(num_layers)]

    def add_parity_bound(self):
        """Adds parity bound."""
        for i in range(self.num_layers // 2):
            low_layer = self.layers[i]
            high_layer = self.layers[self.num_layers - i - 1]
            self.lp.add_constraint(
                [("clique_parity", -1), ((self.v, high_layer), 1), ((self.v, low_layer), -1)],
                ">=",
                self.basis.xor_upper_bound(2))

    def add_counting_bounds(self):
        """Adds counting bounds, for a given number of gates.

        This adds a counting lower bound, for any given number of gates.
        """
        # FIXME


    def add_zeroing_bound(self):
        """Adds bound from zeroing out one vertex."""
        for v in range(self.k, self.n):
            # "at least one NAND gate was hit"
            self.lp.add_constraint(
                [(("V", v+1), 1), (("V", v), -1)],
                ">=",
                1
            )

    def get_all_bounds(self):
        """Gets bounds for each possible number of cliques.

        This is the bounds for each possible number of cliques,
        in the scenario that the number of gates for
        CLIQUE is minimized.
        """
        # solve, minimizing number of gates in "highest layer"
        r = self.lp.solve(("L", self.layers[-1][0]))
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
    parser.add_argument("--max-gates-search",
        help="Search for smallest number of gates which doesn't crash",
        action="store_true")
    parser.add_argument("--result-file",
        help="Write result to indicated file (rather than stdout)")
    return parser.parse_args()





if __name__ == '__main__':
    args = parse_args()
    n = args.n
    k = args.k
    max_gates = args.max_gates

    gate_range = range(args.max_gates, 10000) if args.max_gates_search else [args.max_gates]

    for max_gates in gate_range:
        try:
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
            print(f"***** wrote bound with max_gates={max_gates}")
            sys.exit(0)
        except TypeError:
            print(f"***** failed with max_gates={max_gates}; retrying")

