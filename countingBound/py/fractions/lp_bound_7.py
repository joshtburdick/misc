#!/usr/bin/env python3
# Attempt at LP, using the bouncing walk and exact arithmetic.
# (Note that this doesn't use the LP relaxation of the IP.)

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
import scip_helper

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

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        This will have two groups of variables:
        - tuples of the form ("E", c) where:
            - "c" is the number of cliques, with 0 <= c <= {n choose k}
            - Each variable will be the expected number of gates in
            the sets with exactly c cliques.
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
        # wrapper for LP solver
        self.lp = scip_helper.SCIP_Helper(
            self.expected_num_gates_vars + self.num_gates_dist_vars)

        # number of possible cliques
        self.num_possible_cliques = comb(n, k)
        # number of cliques which include an arbitrary edge: this
        # is both the maximum removed by zeroing, and the maximum added
        self.max_cliques_with_edge = comb(n-2, k-2)

        self.z = self.edge_zeroing_transition_matrix()
        self.a = self.edge_adding_transition_matrix()

        self.basis = gate_basis.UnboundedFanInNandBasis()
        # self.basis = gate_basis.TwoInputNandBasis()

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
                print(f"{i} {j} {z[(i,j)]}")
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
                a[(i,j)] = comb(self.max_cliques_with_edge, j-i) / num_sets
        return a

    def step_probability(self, i, k):
        """Computes probability of stepping from i to k cliques.

        This is basically matrix multiplication, except that we're
        using sparse matrices of fractions, and so using numpy might
        lose precision.
        """
        # loop through possible number of cliques remaining,
        # after zeroing out an edge
        a = fractions.Fraction(0,1)
        for j in range(max(0, i-self.max_cliques_with_edge)+1):
            try:
                a += self.z[(i,j)] * self.a[(j,k)]
            except:
                print(f"couldn't compute step {i} -> {j} -> {k}")
                # pass
        return a

    def add_counting_bound(self):
        """Adds counting bounds, for a given number of gates.

        This is a lower bound on the number of functions with
        some number of gates.
        """
        # lower bound on expected number of gates, for _all_ the functions
        expected_gates_lower_bound = int(self.basis.expected_gates(
            comb(self.n, 2), self.num_possible_cliques))
        # add lower bound on all of these (weighted by number
        # of functions at each level)
        self.lp.add_constraint(
            [(("E", c), comb(self.num_possible_cliques, c))
                for c in range(self.num_possible_cliques+1)],
            ">=", 2**self.num_possible_cliques * expected_gates_lower_bound)

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
            A = [(("E", j), self.step_probability(i, j))
                for j in range(max(0, i-self.max_cliques_with_edge),
                    min(i+self.max_cliques_with_edge, N)+1)]
            A += [(("E", i), fractions.Fraction(-1, 1))]
            # note that for these, we ignore the fact that zeroing out an edge
            # usually removes one gate
            if lower_bound:
                # note that the lower bound for A_{i+1} is actually _less_ than A_i !
                # ummm... it is still a lower bound. uh... it's not clear that this will help.
                self.lp.add_constraint(A, ">=",
                    - fractions.Fraction(i, self.num_possible_cliques) * self.max_cliques_with_edge)
            if upper_bound:
                self.lp.add_constraint(A, "<=",
                    fractions.Fraction(self.max_cliques_with_edge, 2))

    def get_all_bounds(self):
        """Gets bounds, minimizing each variable individually."""
        def get_bound_for_variable(x):
            (num_vertices, num_cliques) = x
            # print(f"solving: n={num_vertices} num_cliques={num_cliques}")
            r = self.lp.solve(x)
            # XXX hack to deal with if this fails
            # (although 0. is legit a lower bound)
            if not r:
                r = 0.
            return (num_vertices, num_cliques, np.round(r, 2))
        # we only get bounds for "expected number of gates"
        bounds = [get_bound_for_variable(x)
            for x in self.expected_num_gates_vars]
        return pandas.DataFrame(bounds,
            columns = ['Num. vertices', 'Num. cliques', 'Min. gates'])

def get_bounds(n, k, constraints_label,
        use_counting_bound, use_lower_bound, use_upper_bound):
    """Gets bounds with some set of constraints.

    n, k: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_lower_bound, use_upper_bound: whether to use
        each of these groups of constraints
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, label={constraints_label}]\n')
    bound = LpEdgeZeroing(n, k)
    if use_counting_bound:
        bound.add_counting_bound()
    if use_lower_bound or use_upper_bound:
        bound.add_step_constraints(use_lower_bound, use_upper_bound)
    b = bound.get_all_bounds()
    b['Constraints'] = constraints_label
    return b.iloc[:,[3,0,1,2]]

if __name__ == '__main__':
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    bounds = pandas.concat([
        get_bounds(n, k, 'Counting', True, False, False),
        get_bounds(n, k, 'Step', False, True, True),
        get_bounds(n, k, 'Counting and step', True, True, True),
    ])
    bounds.to_csv(sys.stdout, index=False)

