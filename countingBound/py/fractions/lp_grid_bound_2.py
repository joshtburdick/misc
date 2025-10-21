#!/usr/bin/env python3
# LP bound, using a "grid" of expected values.

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

class LpEdgeZeroing:
    """Attempt at bound by zeroing out edges.
    """

    def __init__(self, n, k, max_gates):
        """Constructor gets graph info, and sets up variable names.

        This will have two groups of variables, for expected number
        of gates:
        - tuples of the form (i, j), where
            - "i" is the number of cliques _not_ hit by edge "e"
            - "j" is the number of cliques hit by edge "e"
        - tuples of the form ("E", c) where:
            - "c" is the number of cliques, with 0 <= c <= {n choose k}
              (this is the average for a given line of i+j, in the
              above variables)
        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')

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
        # FIXME: make this an option?
        self.lp = pulp_helper.PuLP_Helper(
             self.expected_num_gates_vars + self.num_gates_dist_vars)

        # FIXME make this an option?
        self.basis = gate_basis.UnboundedFanInNandBasis()
        # self.basis = gate_basis.TwoInputNandBasis()

        # number of possible cliques
        self.N = comb(n, k)
        # maximum number of cliques which could be "hit" by an arbitrary edge
        self.max_cliques_hit = comb(n-2, k-2)
        # maximum number of cliques which could be "missed" by that edge
        self.max_cliques_missed = self.N - self.max_cliques_hit

        # number of functions in each set
        self.num_functions = np.zeros(
            (self.max_cliques_missed, self.max_cliques_hit), dtype=int)
        for i in range(self.max_cliques_missed+1):
            for j in range(self.max_cliques_hit+1):
                self.num_functions[i,j] = (
                    comb(self.max_cliques_missed, i)
                    * comb(self.max_cliques_hit, j))
        # the number of functions should add up to this
        assert(self.num_functions.sum() == 2 ** self.N)

    def add_level_constraints(self):
        """Adds average for each total number of cliques."""
        # k is the total number of cliques
        for k in range(self.N+1):
            A = []
            total_functions = 0
            # j is the number of cliques "hit"
            for j in range(self.max_cliques_hit+1):
                # i is the number of cliques "missed"
                i = k - j
                if not (0 <= i <= self.max_cliques_missed):
                    continue
                A += ((i, j), self.num_functions[i, j])
                total_functions += self.num_functions[i, j]
            self.lp.add_constraint(A + [(("E", k), -num_functions)],
                '=', 0)

    def add_step_constraints(self, lower_bound=True, upper_bound=True):
        """Adds 'step' constraints, for tweaking one edge.

        This bounds the number of gates, after one 'step' of zeroing
        an edge.
        """ 
        for i in range(self.num_cliques_missed+1):
            for j in range(1, self.num_cliques_hit+1):
                if lower_bound:
                    # zonking the edge removes at least one gate
                    self.lp.add_constraint(
                        [((i, j), 1), ((i, 0), -1)],
                        ">=", 1)
                if upper_bound:
                    # number of additional gates we'd need, to also
                    # detect the cliques hitting the edge
                    num_additional_gates = (
                        self.basis.or_of_and_upper_bound(comb(k, 2), j)
                        + self.basis.or_upper_bound())
                    self.lp.add_constraint(
                        [((i, j), 1), ((i, 0), -1)],
                        "<=", num_additional_gates)

    def add_counting_bound(self):
        """Adds counting bound.

        This is a weighted average of all of the sets.
        """
        A = []
        for i in range(self.num_cliques_missed+1):
            for j in range(self.num_cliques_hit+1):
                A += [((i,j), self.num_functions([i,j])
        self.lp.add_constraint(A, ">=",
            self.basis.expected_gates(comb(self.n, 2), self.N) /
            (2 ** self.N))

    def add_no_cliques_constraint(self):
        """Adds trivial constraint, on finding no cliques."""
        # ... that is, finding zero cliques requires zero NAND gates
        self.lp.add_constraint([((0, 0), 1)], "=", 0)

    # FIXME add combining bound
    def add_combining_bound(self):
        # starting with a circuit which detects i cliques
        for i in range(1, self.N):
            # ... and a circuit which detects no more than i cliques
            for j in range(1, i+1):
                # we can combine them to detect more cliques
                for k in range(i+1, min(i+j, self.N)+1):
                self.lp.add_constraint(
                    [(("E",k),1), (("E",i),-1), (("E",j),-1)]
                    "<=", self.basis.or_bound())

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
        use_counting_bound, use_lower_bound, use_upper_bound,
        use_no_cliques_bound):
    """Gets bounds with some set of constraints.

    n, k, max_gates: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_lower_bound, use_upper_bound, use_no_cliques_bound:
        whether to use each of these groups of constraints
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, max_gates={max_gates}, label={constraints_label}]\n')
    bound = LpEdgeZeroing(n, k, max_gates)
    bound.add_level_constraints()
    if use_counting_bound:
        bound.add_counting_bound()
    if use_lower_bound or use_upper_bound:
        bound.add_step_constraints(use_lower_bound, use_upper_bound)
    if use_no_cliques_bound:
        bound.add_no_cliques_constraint()
        # XXX this is somewhat misnamed
        bound.add_simple_step_constraint()
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
        get_bounds(n, k, max_gates, 'Step', False, True, True, False),
        get_bounds(n, k, max_gates, 'Counting and step', True, True, True, False),
        get_bounds(n, k, max_gates, 'Counting, step, and no cliques',
            True, True, True, True)
    ])
    if args.result_file:
        with open(args.result_file, "wt") as f:
            bounds.to_csv(f, index=False)
    else:
        bounds.to_csv(sys.stdout, index=False)

