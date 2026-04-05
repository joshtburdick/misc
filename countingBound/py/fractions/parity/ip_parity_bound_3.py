#!/usr/bin/env python3
# IP bound for clique parity.

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

sys.path.append("..")   # XXX

import gate_basis
import hypergraph_counter
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

class LpParity:
    """Attempt at bound for the clique parity problem.
    """

    def __init__(self, n, k, max_gates):
        """Constructor gets graph info, and sets up variable names.

        This will have several groups of variables:
        - tuples of the form ("X", c) where:
            - "c" is the number of cliques, with 0 <= c <= {n choose k}
            - Each variable will be the expected number of gates in
            the sets with `c` cliques.
        - tuples of the form ("A", c, v) where:
            - "c" is the number of cliques, with 0 <= c <= {n choose k}
            - "v" is the number of vertices, with k <= v <= n
            - Each variable will be the expected number of gates in
            the sets with `c` cliques and _exactly_ `v` vertices.
        - tuples of the form ("G", v, g), where
            - "v" is the number of vertices, with k <= v <= n
            - "g" is the number of gates, with 1 <= g <= max_gates
            - Each variable will be the number of sets with _up to_ `v` vertices,
              having _exactly_ `g` gates.
        - tuples of the form ("U", v), which are the expected number of gates
            in functions with _up to_ `v` vertices
        - tuples of the form ("V", v), which, similarly, are the expected
            number of gates, but in functions with _exactly_ `v` vertices
        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        max_gates: maximum number of gates to consider (setting this too
            small will result in an infeasible LP)
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        self.max_gates = max_gates

        self.vars = []
        for c in range(0, comb(n, k) + 1):
            # variables for expected number of gates, for each number of cliques
            self.vars += [("X", c)]
        for v in range(k, n+1):
            # averages by number of vertices and cliques
            for c in range(0, comb(n, k) + 1):
                self.vars += [("A", v, c)]
            # averages by number of vertices
            # ("U" is "up to", "V" is "exactly", v vertices)
            self.vars += [("U", v), ("V", v)]
            # counts of numbers of functions with some number of gates
            # (with up to some number of vertices)
            for g in range(1, max_gates+1):
                self.vars += [("G", v, g)]

        # wrapper for LP solver
        self.lp = pulp_helper.PuLP_Helper(self.vars)
 
        # number of possible cliques
        self.num_possible_cliques = comb(n, k)

        self.basis = gate_basis.TwoInputNandBasis()

        # counts of hypergraphs with _exactly_ some number of vertices
        self.hc = hypergraph_counter.HypergraphCounter(self.n, self.k)
        self.hypergraph_counts = self.hc.count_hypergraphs_exact_vertices()

        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None

    def add_level_constraints(self):
        """Adds constraints on functions at some "level".

        By "level", we mean "number of vertices".

        The constraints (both of which are equalities) are:
        - on the total number of functions at that "level", and
        - connecting the counts with that "level"'s expected gate count
        """
        # loop through number of vertices; note that we're counting functions
        # with _exactly v_ vertices (as opposed to _up to v_ vertices).
        for v in range(self.k, self.n+1):
            num_functions = self.hypergraph_counts[v].sum()
            # add constraint that these sum to the number of functions
            self.lp.add_constraint(
                [(("G", v, i), 1) for i in range(1, self.max_gates+1)],
                '=', num_functions)
            # add constraint defining expected number of gates
            A = [(("G", v, i), i)
                for i in range(1, self.max_gates+1)]
            self.lp.add_constraint(A + [(("V", v), -num_functions)],
                '=', 0)

    def add_cumulative_constraints(self):
        """Adds constraints that connect the "up to v" and "exactly v" variables.
        """
        num_hypergraphs_by_num_vertices = {
            v: counts.sum()
            for v, counts in self.hypergraph_counts.items()
        }
        # This will be the total number of hypergraphs with up to and including
        # v vertices.
        total_hypergraphs = 0
        for v in range(self.k, self.n+1):
            total_hypergraphs += num_hypergraphs_by_num_vertices[v]
            # The average number of gates in functions with up to v vertices
            # is the weighted average of the average number of gates in functions
            # with exactly i vertices, for k <= i <= v.
            self.lp.add_constraint(
                [(("V", i), num_hypergraphs_by_num_vertices[i]) for i in range(self.k, v+1)]
                + [(("U", v), -total_hypergraphs)],
                '=', 0)

    def add_marginal_constraints(self):
        """Adds constraints on some 'marginal' variables.

        """
        # "V", the number of gates needed to detect functions with exactly some
        # number of vertices, is a weighted average of the numbers in
        # "A", which is the number of gates needed, grouped by _exact_ number
        # of vertices and number of cliques.
        for v in range(self.k, self.n+1):
            n_hypergraphs = self.hypergraph_counts[v]
            self.lp.add_constraint(
                [(("A", v, c), n_hypergraphs[c]) for c in range(n_hypergraphs.shape[0])]
                + [(("V", v), -n_hypergraphs.sum())],
                '=', 0)
        # "X", the number of gates needed to detect functions with exactly some
        # number of cliques, similarly is a weighted average of numbers in "A".
        for c in range(self.num_possible_cliques+1):
            n_functions = {
                v: self.hypergraph_counts[v][c]
                for v in range(self.k, self.n+1)
                if c < self.hypergraph_counts[v].shape[0]
            }
            # ??? is this right?
            self.lp.add_constraint(
                [(("A", v, c), f) for v, f in n_functions.items()]
                + [(("X", c), -sum(n_functions.values()))],
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
                [(("G", v, g), 1)
                    for v in range(self.k, self.n+1)],
                '<=', num_functions[g])

    def add_zeroing_constraint(self):
        """Constraint based on zeroing out one vertex.

        Note that we use the "up to v cliques" variables here. This is
        because after we zero out one vertex, we don't know how many
        vertices remain.
        """
        # number of vertices in the larger graph
        for v in range(self.k+1, self.n+1):
            # number of cliques "hit" by zeroing out edges connected to one vertex
            num_hit_cliques = comb(v-1, self.k-1)
            # probability that we hit at least one clique
            prob_hit_at_least_one_clique = 1 - 2 ** (-num_hit_cliques)
            # Note that if we hit one clique (as we usually do), we _only_ can
            # say that one NAND gate was zonked.
            self.lp.add_constraint(
                [(("U", v), 1), (("U", v-1), -1)],
                ">=", prob_hit_at_least_one_clique)

    def add_one_clique_constraint(self):
        """Constraint on the number of gates to detect one clique."""
        # Constraint for parity of zero cliques (??? is this right?)
        self.lp.add_constraint([(("E", 0), 1)], "=", 1)
        # FIXME currently this is hard-coded for unbounded fan-in
        self.lp.add_constraint([(("E", 1), 1)], "=", 2)

    def add_smoothing_constraint(self):
        """Constraint on computing parity with one clique added or removed."""
        # Number of gates to detect one clique, and XOR it with the rest.
        # For unbounded fan-in NAND gates, this is 5: 2 for the clique,
        # and 3 for the XOR (as it can reuse one one-input NAND gate).
        delta_gates = 5
        for i in range(self.num_possible_cliques):
            self.lp.add_constraint(
                [(("E", i+1), 1), (("E", i), -1)],
                "<=", delta_gates)
            self.lp.add_constraint(
                [(("E", i+1), 1), (("E", i), -1)],
                ">=", -delta_gates)

    def add_all_cliques_constraint(self):
        """Adds constraint on computing parity of all cliques.

        If we have a way to compute parity of all N cliques, we can
        use it, plus an XOR, to convert a circuit which computes
        parity of i cliques, into one which computes parity of N-i cliques.

        E_N >= E_{N-i} - E_i + 4
        E_{N-i} - E_i - E_n <= -4
        """
        N = self.num_possible_cliques
        for i in range(1, N // 2):
            self.lp.add_constraint(
                [(("E", N), -1), (("E", N-i), 1), (("E", i), -1)],
                "<=", -4)

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

