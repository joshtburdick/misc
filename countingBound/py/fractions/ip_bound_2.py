#!/usr/bin/env python3
# Revised IP bound.
# This is intended to allow using an IP, or the LP relaxation.

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
import exact_simplex_helper

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
        # FIXME: make this an option?
        # self.lp = scip_helper.SCIP_Helper(
        #     self.expected_num_gates_vars + self.num_gates_dist_vars)
        # self.lp = flexible_lp_helper.Flexible_LP_Helper(
        #      self.expected_num_gates_vars + self.num_gates_dist_vars)
        self.lp = pulp_helper.PuLP_Helper(
             self.expected_num_gates_vars + self.num_gates_dist_vars)
        # self.lp = exact_simplex_helper.ExactSimplexHelper(
        #     self.expected_num_gates_vars + self.num_gates_dist_vars)
        # number of possible cliques
        self.num_possible_cliques = comb(n, k)
        # number of cliques which include an arbitrary edge: this
        # is both the maximum removed by zeroing, and the maximum added
        self.max_cliques_with_edge = comb(n-2, k-2)

        self.z = self.edge_zeroing_transition_matrix()
        self.a = self.edge_adding_transition_matrix()

        self.basis = gate_basis.UnboundedFanInNandBasis()
        # self.basis = gate_basis.TwoInputNandBasis()
        self.rng = np.random.default_rng()

        # for debugging: directory in which to save LP problem files
        self.lp_save_dir = None

    def random_eps(self):
        """Generates a small number, to perturb constraints a bit.

        Hoping this may avoid some crashes...
        """
        # return fractions.Fraction(
        #     self.rng.integers(1, 100),
        #     100000000)
        # currently disabled
        return 0

    def edge_zeroing_transition_matrix(self):
        """The zeroing transition matrix.

        This is for zeroing out one edge.
        This follows the convention for Markov chains in Wikipedia:

        - each of its rows sums to 1, so it would be
          "right stochastic", except that it's rectangular
        - if x is a row vector of probabilities, then xZ is the
          row vector of probabilities after one step of zeroing
        """
        N = comb(self.n, self.k)
        # we make the matrix rectangular
        z = dict()
        for i in range(N+1):
            for j in range(max(0, i-self.max_cliques_with_edge),
                    min(i+1, N+1-self.max_cliques_with_edge)):
                z[(i,j)] = hyperg_frac(N, i, N-self.max_cliques_with_edge, j)
                # print(f"{i} {j} {z[(i,j)]}")
        return z

    def edge_adding_transition_matrix(self):
        """The edge-adding transition matrix.

        This follows the convention for Markov chains in Wikipedia:

        - each of its rows sums to 1, so it would be
          "right stochastic", except that it's rectangular
        - if x is a row vector of probabilities, then xZ is the
          row vector of probabilities after one step of zeroing
        """
        N = comb(self.n, self.k)
        # number of possible sets of cliques which could be added
        num_sets = 2 ** self.max_cliques_with_edge 
        # we make the matrix rectangular
        a = dict()
        for i in range(N+1-self.max_cliques_with_edge):
            for j in range(i, min(i+1+self.max_cliques_with_edge, N+1)):
                # ??? is this right?
                a[(i,j)] = fractions.Fraction(comb(self.max_cliques_with_edge, j-i), num_sets)
        return a

    def step_probability(self, i, k):
        """Computes probability of stepping from i to k cliques.

        This is basically matrix multiplication, except that we're
        using sparse matrices of fractions, and so using numpy might
        lose precision.
        """
        # loop through possible number of cliques remaining,
        # after zeroing out an edge
        s = fractions.Fraction(0)
        # FIXME only looping over a subset of these isn't working
        # (but would be faster)
        # was: for j in range(max(0, i-self.max_cliques_with_edge)+1, i+1):
        for j in range(max(0, i-self.max_cliques_with_edge),
            min(i+self.max_cliques_with_edge, self.num_possible_cliques)+1):
        # OK, this is definitely wrong...
        # for j in range(self.max_cliques_with_edge+1):
            try:
                s += self.z[(i,j)] * self.a[(j,k)]
            except:
                # We treat missing values here as 0's, relying on
                # the sparsity of Z and A here. (Thus, if one of these
                # is missing, it shouldn't be a problem.)
                # print(f"couldn't compute step {i} -> {j} -> {k}")
                pass
        return s

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
                '=', num_functions + self.random_eps())
            # add constraint defining expected number of gates
            A = [((c, i), i)
                for i in range(1, self.max_gates+1)]
            self.lp.add_constraint(A + [(("E", c), -num_functions)],
                '=', 0 + self.random_eps())

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
                '<=', num_functions[g] + self.random_eps())

    def get_step_bounds(self):
        """Gets the vectors of change in number of gates at each step.

        `num_added[i]` is the expected number of cliques which _were just added_ to obtain
        a set of `i` cliques, in the "up" step.

        Somewhat symmetrically, `num_hit[i]` is the number of cliques which were "hit" in the
        "down" step, obtaining a set of `i` cliques. Since we're counting this "after" the bounce,
        the numbers are "blurred" by the transition matrix.
        """
        num_added = [fractions.Fraction(self.num_cliques_with_edge * i, self.num_possible_cliques)
            for i in range(self.num_possible_cliques + 1)]
        num_hit = [-sum([self.step_probability(i, j) * num_added(j) for j in range(self.num_possible_cliques+1)])
            for i in range(self.num_possible_cliques+1)]
        return (num_hit, num_added)

    def add_step_constraints(self, lower_bound=True, upper_bound=True):
        """Adds 'step' constraints, for tweaking one edge.

        This bounds the number of gates, after one 'step' of zeroing out
        an edge, and then adding in cliques incident to it.
        That is, it bounds |C(A_{i+1})|, relative to |C(A_i)|.
        lower_bound: if True, then include the lower bound.
        upper_bound: if True, then include the upper bound.
        """
        N = self.num_possible_cliques
        for i in range(N+1):
#            A = [(("E", j), self.step_probability(i, j))
#                for j in range(max(0, i-self.max_cliques_with_edge),
#                    min(i+self.max_cliques_with_edge, N)+1)]
# XXX for now, just including the entire row of LU
            A = [(("E", j), self.step_probability(i, j))
                for j in range(self.num_possible_cliques+1)]
            A += [(("E", i), fractions.Fraction(-1))]
            # pdb.set_trace()
            # print(A)
            # note that for these, we ignore the fact that zeroing out an edge
            # usually removes one gate
            if lower_bound:
                # note that the lower bound for A_{i+1} is actually _less_ than A_i !
                # ummm... it is still a lower bound. uh... it's not clear that this will help.
                self.lp.add_constraint(A, ">=",
                    - (fractions.Fraction(i, self.num_possible_cliques) * self.max_cliques_with_edge) + fractions.Fraction(0))
            if upper_bound:
                self.lp.add_constraint(A, "<=",
                    fractions.Fraction(self.max_cliques_with_edge, 2)+fractions.Fraction(0))

    def add_no_cliques_constraint(self):
        """Adds trivial constraint, on finding no cliques."""
        # ... that is, finding zero cliques requires one NAND gate
        self.lp.add_constraint([(("E", 0), 1)], "=", 1)

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

def write_transition_matrix_check(n, k):
    """Prints a quick check of the transition matrices."""
    bound = LpEdgeZeroing(n,k,2)
    def dict_to_matrix(x):
        """Converts dict-of-keys to matrix (for printing)."""
        m = max([i for (i,_) in x.keys()]) + 1
        n = max([j for (_,j) in x.keys()]) + 1
        a = scipy.sparse.dok_matrix((m, n), float)
        for ((i,j),x1) in x.items():
            a[i,j] = float(x1)
        return a.todense()
    np.set_printoptions(suppress=True, linewidth=100000)
    with open(f"transition_matrix_check_{n}_{k}.txt", "wt") as f:
        z = dict_to_matrix(bound.z)
        a = dict_to_matrix(bound.a)
        f.write(f"Z =\n{np.round(z, 3)}\n\n")
        f.write(f"A =\n{np.round(a, 3)}\n\n")
        f.write(f"ZA =\n{np.round(z @ a, 3)}\n\n")

        # comparing this with version computed with
        # "written-out" matrix multiplication
        za1 = scipy.sparse.dok_matrix(
            (bound.num_possible_cliques+1, bound.num_possible_cliques+1),
            float)
        for i in range(bound.num_possible_cliques+1):
            for k in range(bound.num_possible_cliques+1):
                za1[i,k] = bound.step_probability(i, k)
        za1 = dict_to_matrix(za1)
        f.write(f"ZA (alternate take) =\n{np.round(za1, 3)}\n\n")

        N = bound.num_possible_cliques
        x = np.array([comb(N,i) for i in range(N+1)])
        f.write(f"binomial coefficients: x = \n{x}\n\n")
        f.write(f"x * Z * A =\n{ x @ z @ a }\n\n")

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
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    n = args.n
    k = args.k
    max_gates = args.max_gates

    # possibly just dump the transition matrix
    if args.write_transition_matrices:
        write_transition_matrix_check(n, k)
        sys.exit(0)

    # output_file = sys.argv[3]
    bounds = pandas.concat([
        get_bounds(n, k, max_gates, 'Counting', True, False, False, False),
        get_bounds(n, k, max_gates, 'Step', False, True, True, False),
        get_bounds(n, k, max_gates, 'Counting and step', True, True, True, False),
        get_bounds(n, k, max_gates, 'Counting, step, and no cliques',
            True, True, True, True)
    ])
    bounds.to_csv(sys.stdout, index=False)

