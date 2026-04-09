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
sys.path.append("../../")   # XXX

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
        (FIXME rename these somewhat?)
        - tuples of the form ("X", c) where:
            - "c" is the number of cliques, with 0 <= c <= {n choose k}
            - Each variable will be the expected number of gates in
            the sets with `c` cliques.
        - tuples of the form ("A", c, v) where:
            - "c" is the number of cliques, with 0 <= c <= {v choose k}
            - "v" is the number of vertices, with k <= v <= n
            - Each variable will be the expected number of gates in
            the sets with `c` cliques and _exactly_ `v` vertices.
        - tuples of the form ("G", c, v, g), where
            - "c" is the number of cliques, with 0 <= c <= {v choose k}
            - "v" is the number of vertices, with k <= v <= n
            - "g" is the number of gates, with 1 <= g <= max_gates
            - Each variable will be the number of sets with `c` cliques,
              and _exactly_ `v` vertices, having _exactly_ `g` gates.
        - tuples of the form ("U", v), which are the expected number of gates
            in functions with _up to_ `v` vertices
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
            # Variables for expected number of gates, for each number of cliques.
            self.vars += [("X", c)]
        for v in range(k, n+1):
            # Averages by number of vertices and cliques. (When the number of vertices is large,
            # the first few of these variables will be zero; we add them anyway.)
            for c in range(0, comb(v, k) + 1):
                self.vars += [("A", c, v)]
                # Counts of numbers of functions with some number of gates
                # (with exactly `v` vertices and `c` cliques).
                for g in range(1, max_gates+1):
                    self.vars += [("G", c, v, g)]
            #   Averages by number of vertices.
            self.vars += [("U", v)]

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
        FIXME rename to "add_set_constraints"? or something like that?

        The constraint (which is an equality) is on the number of functions
        with some number of vertices and cliques.
        """
        # loop through number of vertices; note that we're counting functions
        # with _exactly v_ vertices (as opposed to _up to v_ vertices).
        for v in range(self.k, self.n+1):
            for c in range(0, comb(v, k) + 1):
                num_functions = self.hypergraph_counts[v][c]
                # add constraint that these sum to the number of functions
                self.lp.add_constraint(
                    [(("G", c, v, g), 1) for g in range(1, self.max_gates+1)],
                    '=', num_functions)

    def add_marginal_constraints(self):
        """Adds constraints on some 'marginal' variables.

        """
        # "U", the number of gates needed to detect functions with up to some
        # number of vertices, is a weighted average of the numbers in
        # "A", which is the number of gates needed, grouped by _exact_ number
        # of vertices and number of cliques.
        total_hypergraphs = 0
        A = []
        for v in range(self.k, self.n+1):
            # We add in the number of hypergraphs with exactly v vertices,
            # and any number of cliques (since `U` is for "up to v vertices").
            n_hypergraphs = self.hypergraph_counts[v]
            total_hypergraphs += n_hypergraphs.sum()
            A += [(("A", c, v), n_hypergraphs[c]) for c in range(n_hypergraphs.shape[0])]
            self.lp.add_constraint(A + [(("U", v), -total_hypergraphs)],
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
                [(("A", c, v), f) for v, f in n_functions.items()]
                + [(("X", c), -sum(n_functions.values()))],
                '=', 0)

    def add_counting_bound(self):
        """Adds counting bounds, for a given number of gates.

        This is an upper bound on the number of functions with
        some number of gates.
        """
        # number of possible functions for each possible number of gates
        # (with number of inputs based on number of vertices)
        num_functions = self.basis.num_functions(comb(self.n, 2), self.max_gates+1)
        # The entries in `G` are essentially the counts which go in to computing
        # the averages `A`. Therefore, find variables which are named `A`.
        A_vars = [v for v in self.vars if v[0] == "A"]
        # for each number of gates, upper-bound "total number of functions with this many gates"
        for g in range(1, self.max_gates+1):
            self.lp.add_constraint(
                [(("G", c, v, g), 1) for _, c, v in A_vars],
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
        self.lp.add_constraint([(("X", 0), 1)], "=", 1)
        # FIXME currently this is hard-coded for unbounded fan-in
        self.lp.add_constraint([(("X", 1), 1)], "=", 2)

    def add_smoothing_constraint(self):
        """Constraint on computing parity with one clique added or removed."""
        # Number of gates to detect one clique, and XOR it with the rest.
        # For unbounded fan-in NAND gates, this is 5: 2 for the clique,
        # and 3 for the XOR (as it can reuse one one-input NAND gate).
        delta_gates = 5
        for i in range(self.num_possible_cliques):
            # E_{i+1} <= E_i + 5, so
            # E_{i+1} - E_i <= 5
            self.lp.add_constraint(
                [(("X", i+1), 1), (("X", i), -1)],
                "<=", delta_gates)
            # E_{i+1} >= E_i - 5, so
            # E_{i+1} - E_i >= -5
            self.lp.add_constraint(
                [(("X", i+1), 1), (("X", i), -1)],
                ">=", -delta_gates)

    def add_all_cliques_constraint(self):
        """Adds constraint on computing parity of all cliques.

        If we have a way to compute parity of all N cliques, we can
        use it, plus an XOR, to convert a circuit which computes
        parity of i cliques, into one which computes parity of N-i cliques.

        (It's not clear that this is useful.)

        X_{n-i} <= X_n + X_i + 4
        X_{n-i} - X_n - X_i <= 4
        """
        N = self.num_possible_cliques
        for i in range(1, N // 2):
            self.lp.add_constraint(
                [(("X", N-i), 1), (("X", N), -1), (("X", i), -1)],
                "<=", 4)

    def get_all_bounds(self):
        """Gets bounds for each possible number of cliques.

        This is the bounds for each possible number of cliques,
        in the scenario that the number of gates for
        CLIQUE is minimized.
        """
        # solve, minimizing number of gates for CLIQUE
        r = self.lp.solve(("X", self.num_possible_cliques))
        if not r:
            return None
        # for now, we only get bounds for "expected number of gates"
        # for each number of cliques
        n_cliques = range(self.num_possible_cliques+1)
        bounds = [r[("X", num_cliques)]
            for num_cliques in range(self.num_possible_cliques+1)]
        return pandas.DataFrame({
                'Num. vertices': self.n,
                'Num. cliques': n_cliques,
                'Min. gates': bounds})

def get_bounds(n, k, max_gates, constraints_label,
        use_counting_bound, use_zeroing_bound, use_smoothing_bound):
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
    bound.add_marginal_constraints()
    if use_counting_bound:
        bound.add_counting_bound()
    if use_zeroing_bound:
        bound.add_zeroing_constraint()
    if use_smoothing_bound:
        bound.add_one_clique_constraint()
        bound.add_smoothing_constraint()
        bound.add_all_cliques_constraint()
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
                get_bounds(n, k, max_gates, 'Counting, zeroing', True, True, False),
                get_bounds(n, k, max_gates, 'Counting, smoothing', True, False, True),
                get_bounds(n, k, max_gates, 'Counting, zeroing, smoothing',
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

