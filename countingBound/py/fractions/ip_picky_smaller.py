#!/usr/bin/env python3
# Revised IP bound, using the "picky" version of BUGGYCLIQUE, but smaller.

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
	Here, we consider PICKYCLIQUE to be the cases with at least one YES clique,
	*and* at least one NO clique. (The special case with no NO cliques is
	the previous BUGGYCLIQUE bound.
	"""
    def __init__(self, n, k, max_gates):
        """Constructor gets graph info, and sets up variable names.

        This will have two groups of variables:
        - tuples of the form `(function_type, num_cliques, "E")` where:
			- `function_type` is either "buggy" or "picky"
			- `num_cliques` is the *total* number of YES and NO cliques
			  (for "buggy", this starts at 1; for "picky", this starts at 2)
			- `"E"` means "the expected number of gates"; and
        - tuples of the form `(function_type, num_cliques, num_gates)` where:
			- `function_type` and `num_cliques` are as before
			- `num_gates` (ranging from 1 to the user-defined max. number of gates)	
			  is the number of gates
		Constructor args:

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

        for function_type in ["buggy", "picky"]:
            min_cliques = 1 if function_type=="buggy" else 2
            for i in range(min_cliques, self.num_possible_cliques + 1):
                # expected number of gates, for each number of cliques
                self.expected_num_gates_vars += [(function_type, i, "E")]
                for g in range(1, max_gates+1):
                    # counts of numbers of functions with some number of gates
                    self.num_gates_dist_vars += [(function_type, i, g)]

        # wrapper for LP solver
        # FIXME: make this an option?
        self.lp = pulp_helper.PuLP_Helper(
             self.expected_num_gates_vars + self.num_gates_dist_vars)

        # FIXME make this an option?
        # self.basis = gate_basis.UnboundedFanInNandBasis()
        self.basis = gate_basis.TwoInputNandBasis()

        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None

    def add_num_functions_constraints(self):
        """Adds constraints on number of functions at some "level".

        By "level", we mean "number of cliques".

        The constraints (both of which are equalities) are:
        - on the total number of functions at that "level", and
        - connecting the counts with that "level"'s expected gate count
        """
        # adds constraints for one set of functions
        def add_constraints(function_type, num_cliques, num_functions):
            # add constraint defining expected number of gates
            A = [((function_type, num_cliques, g), g)
                for g in range(1, self.max_gates+1)]
            self.lp.add_constraint(A + [((function_type, num_cliques, "E"), -num_functions)],
                '=', 0)
            # add constraint that these sum to the number of functions
            self.lp.add_constraint(
                [((function_type, num_cliques, g), 1) for g in range(1, self.max_gates+1)],
                '=', num_functions)

        for i in range(1, self.num_possible_cliques+1):
            # number of ways to pick i cliques
            add_constraints("buggy", i, comb(self.num_possible_cliques, i))
        for i in range(2, self.num_possible_cliques+1):
            # number of ways to pick i cliques, times the number of ways of labelling them
            # YES or NO (except not all NO)
            add_constraints("picky", i, comb(self.num_possible_cliques, i) * (2**i-1))

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
            A = ([(("buggy", i, g), 1) for i in range(1, self.num_possible_cliques+1)]
                + [(("picky", i, g), 1) for i in range(2, self.num_possible_cliques+1)])
            self.lp.add_constraint(A, '<=', num_possible_functions[g])

    def add_buggy_bound(self):
        """Adds upper bound on computing 'buggy' sets of functions.

        We can implement BUGGYCLIQUE(A+B) by combining circuits:
        PICKYCLIQUE(A,B) OR BUGGYCLIQUE(D), where B <= D <= A+B.

        Here, we average over all the different sets A, weighted by their counts.
        """
        for i in range(2, self.num_possible_cliques + 1):
            # Number of "buggy" functions in this case, which is not all or none
            # of the relevant cliques. We scale everything by this amount to
            # avoid fractions (which sometimes toasts the solver).
            num_buggy = 2**i - 2
            self.lp.add_constraint(
                [(("buggy",i,"E"), num_buggy), (("picky",i,"E"), -num_buggy)]
                + [(("buggy",j,"E"), -comb(i, j)) for j in range(1, i)],
                "<=",
                self.basis.or_upper_bound())

    def add_picky_bound(self):
        """Adds upper bound on computing 'picky' sets of functions.

        We can implement PICKYCLIQUE(A,B) by combining circuits:
        BUGGYCLIQUE(D) AND NOT BUGGYCLIQUE(B), where A <= D <= A+B.
        """
        for i in range(2, self.num_possible_cliques + 1):
            num_buggy = 2**i - 2
            self.lp.add_constraint(
                [(("picky",i,"E"), num_buggy), (("buggy",i,"E"), -num_buggy)]
                + [(("buggy",j,"E"), -comb(i, j)) for j in range(1, i)],
                "<=",
                self.basis.and_upper_bound() + self.basis.not_upper_bound())
                # FIXME double check "AND NOT" cost?

    def add_naive_upper_bound(self):
        """Adds naive upper bound for finding smallish numbers of cliques."""
        # add bound for finding one or more cliques
        for num_cliques in range(1, self.num_possible_cliques+1):
            num_gates = self.basis.or_of_and_upper_bound(
                comb(self.n, 2), num_cliques)
            if num_gates <= self.max_gates:
                self.lp.add_constraint([(("buggy",num_cliques,"E"), 1)],
                    "<=", num_gates)
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
        r = self.lp.solve(("buggy", self.num_possible_cliques, "E"))
        if not r:
            return None
        # XXX for now, just getting counts for BUGGYCLIQUE
        n_cliques = range(1, self.num_possible_cliques+1)
        bounds = [r[("buggy", num_cliques, "E")]
            for num_cliques in n_cliques]
        return pandas.DataFrame({
            'Num. vertices': self.n,
            'Num. cliques': n_cliques,
            'Min. gates': bounds})

def get_bounds(n, k, max_gates, constraints_label,
        use_counting_bound, use_buggy_bound, use_picky_bound, use_upper_bound):
    """Gets bounds with some set of constraints.

    n, k, max_gates: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_buggy_bound, use_picky_bound, use_upper_bound:
        whether to use each of these groups of constraints
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, max_gates={max_gates}, label={constraints_label}]\n')
    bound = LpPicky(n, k, max_gates)
    bound.add_num_functions_constraints()
    if use_counting_bound:
        bound.add_counting_bound()
    if use_buggy_bound:
        bound.add_buggy_bound()
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

