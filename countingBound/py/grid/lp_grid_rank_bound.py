#!/usr/bin/env python3
# LP bound, using a "grid" of expected "rank" of functions.

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

class LpGridRank:
    """Attempt at bound by zeroing out edges.
    """

    def __init__(self, n, k):
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

        # number of possible cliques
        self.N = comb(n, k)
        # maximum number of cliques which could be "hit" by an arbitrary edge
        self.max_cliques_hit = comb(n-2, k-2)
        # maximum number of cliques which could be "missed" by that edge
        self.max_cliques_missed = self.N - self.max_cliques_hit

        # number of functions in each set
        self.num_functions = np.full(
            (self.max_cliques_missed+1, self.max_cliques_hit+1), None)
        for i in range(self.max_cliques_missed+1):
            for j in range(self.max_cliques_hit+1):
                self.num_functions[i,j] = (
                    comb(self.max_cliques_missed, i)
                    * comb(self.max_cliques_hit, j))
        # the number of functions should add up to this
        assert(self.num_functions.sum() == 2 ** self.N)

        # define variables
        # "grid" of variables, for each number of cliques "missed" and "hit"
        grid_variables = []
        for i in range(self.max_cliques_missed+1):
            for j in range(self.max_cliques_hit+1):
                grid_variables += [(i,j)]
        # variables for expected number of gates, for each number of cliques
        num_gates_variables = [("E", c)
            for c in range(0, comb(n, k) + 1)]

        # wrapper for LP solver
        # FIXME: make this an option?
        self.lp = pulp_helper.PuLP_Helper(
            grid_variables + num_gates_variables)

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
                A += [((i, j), self.num_functions[i, j])]
                total_functions += self.num_functions[i, j]
            self.lp.add_constraint(A + [(("E", k), -total_functions)],
                '=', 0)

    def add_counting_bound(self):
        """Adds counting bound.

        This is a weighted average of all of the sets.
        """
        A = []
        # these are all weighted by the number of functions
        for i in range(self.max_cliques_missed+1):
            for j in range(self.max_cliques_hit+1):
                A += [((i,j), self.num_functions[i,j] / (2 ** self.N))]
        # the average rank of all of the functions
        self.lp.add_constraint(A, "=", ((2 ** self.N) - 1) / 2)

    def add_step_constraints(self, lower_bound=True, upper_bound=True):
        """Adds 'step' constraints, for tweaking one edge.

        This bounds the rank, after one 'step' of zeroing an edge.
        """ 
        for i in range(self.max_cliques_missed+1):
            FIXME
            A = [((i,j), comb(  , j))
                for j in range(self.max_cliques_hit+1)]
            self.lp.add_constraint(A,
                ">=",
                (comb(N,i)-1) / 2)

    def get_all_bounds(self):
        """Gets bounds for each possible number of cliques.

        This is the bounds for each possible number of cliques,
        in the scenario that the number of gates for
        CLIQUE is minimized.
        """
        # solve, minimizing number of gates for CLIQUE
        r = self.lp.solve(("E", self.N))
        if not r:
            return None
        # for now, we only get bounds for "expected number of gates"
        # for each number of cliques
        n_cliques = range(self.N+1)
        bounds = [r[("E", num_cliques)]
            for num_cliques in range(self.N+1)]
        return pandas.DataFrame({
                'Num. vertices': self.n,
                'Num. cliques': n_cliques,
                'Min. gates': bounds})

def get_bounds(n, k, constraints_label,
        use_counting_bound, use_lower_bound, use_upper_bound,
        use_combining_bound):
    """Gets bounds with some set of constraints.

    n, k: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_lower_bound, use_upper_bound, use_combining_bound:
        whether to use each of these groups of constraints
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, label={constraints_label}]\n')
    bound = LpEdgeZeroing(n, k)
    bound.add_level_constraints()
    if use_counting_bound:
        bound.add_counting_bound()
        # we include this as well, as it seems pretty basic
        bound.add_no_cliques_constraint()
    if use_lower_bound or use_upper_bound:
        bound.add_step_constraints(use_lower_bound, use_upper_bound)
    if use_combining_bound:
        bound.add_combining_bound()

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

    # output_file = sys.argv[3]
    bounds = pandas.concat([
        get_bounds(n, k, 'Counting', True, False, False, False),
        get_bounds(n, k, 'Counting and step', True, True, True, False),
        get_bounds(n, k, 'Counting, step, and combining',
            True, True, True, True)
    ])
    if args.result_file:
        with open(args.result_file, "wt") as f:
            bounds.to_csv(f, index=False)
    else:
        bounds.to_csv(sys.stdout, index=False)

