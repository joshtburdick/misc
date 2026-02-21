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

def group_by_key(items, key_func):
    """This groups items by a key function, without sorting."""
    items = sorted(items, key=key_func)
    return itertools.groupby(items, key=key_func)

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
        self.layers = self.make_layers(num_layers)

        # Number of functions for each group (num. vertices, layer).
        self.counts_by_group = self.get_counts_by_group()

        # The variables for the LP.
        vars = []
        # Loop through the (num. vertices, layer) pairs.
        for ((v, layer, group_size) in self.counts_by_group.items()):
            # Counts number of gates, for functions in each group.
            vars += [((v, layer, g) for g in range(self.max_gates+1))]
        # Expected number of gates, for each group.
        for v, layer in self.counts_by_group:
            vars += [("E", v, layer)]
        # Averages, by number of vertices, and layer
        for v in range(self.k, self.n+1):
            vars += [("V", v)]
        for l1, l2 in self.layers:
            vars += [("L", l1)]
        # wrapper for LP solver
        self.lp = pulp_helper.PuLP_Helper(vars)
        # basis for gates
        self.basis = gate_basis.TwoInputNandBasis()
        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None

    def make_layers(self, num_layers):
        """Makes the layers for the LP.
        
        We divide the possible number of cliques into layers, and add
        constraints on the average number of gates for each layer.
        We require the layers to be symmetric around N/2.
        (The layers won't, in general, be the same height; but this
        tries to make them as similar as possible.)

        FIXME: the "layers balanced about N/2" may not actually be necessary;
        it would be nice to simplify this.

        num_layers: currently must be odd.
        """
        if num_layers % 2 == 0:
            raise ValueError('num_layers must be odd')
        height = self.max_cliques // num_layers
        endpoints = range(0, self.max_cliques, height)
        endpoints += [self.max_cliques-i for i in reversed(endpoints)]
        return [(endpoints[i], endpoints[i+1]) for i in range(num_layers)]

    def get_counts_by_group(self):
        """Gets number of functions, for each possible (num. vertices, layer).
        
        Returns a dictionary mapping (num. vertices, layer) to number of functions.
        Note that many of the groups will be empty, and so this will not include
        all possible (num. vertices, layer) pairs.
        """
        hypergraph_counts = hypergraph_counts.HypergraphCounts(self.n, self.k)
        counts_by_num_vertices = hypergraph_counts.count_hypergraphs_exact_vertices()
        counts_by_group = {}
        for v in range(self.k, self.n+1):
            for l1, l2 in self.layers:
                if l1 <= len(counts_by_num_vertices[v]):
                    counts_by_group[(v, l1)] = counts_by_num_vertices[v][l1..l2].sum()
        # check that these add up to the total number of functions
        assert sum(counts_by_group.values()) == 2 ** self.num_possible_cliques
        return counts_by_group

    def add_averaging_constraints(self):
        """Adds constraints on the average number of gates for each group."""
        for v, layer in self.counts_by_group:
            A = [((v, layer, g), g) for g in range(self.max_gates+1)]
            # "expected number of gates" = sum(counts * gates) / sum(counts)
            self.lp.add_constraint(
                [(("E", v, layer), -1)] + A,
                "=",
                self.counts_by_group[(v, layer)]
            )
            # add constraints on the number of functions in each group
            A = [((v, layer, g), 1) for g in range(self.max_gates+1)]
            self.lp.add_constraint(
                "=",
                self.counts_by_group[(v, layer)]
            )

    def add_counting_bounds(self):
        """Adds counting bounds, for a given number of gates."""
        num_possible_functions = self.basis.num_functions(
            comb(self.n, 2),
            max_gates
        )
        for g in range(self.max_gates+1):
            self.lp.add_constraint(
                [((v, layer, g), 1) for v, layer in self.counts_by_group],
                "<=",
                num_possible_functions[g]
            )

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
    if use_zeroing_bound:
        bound.add_zeroing_bound()

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
    parser.add_argument("num_layers", type=int,
        help="Number of layers")
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
    num_layers = args.num_layers
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
        # XXX
        except TypeError:
            print(f"***** failed with max_gates={max_gates}; retrying")
