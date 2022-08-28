#!/usr/bin/env python3
# Bounds the rank of functions, in the full lattice, using linear programming.
# FIXME
# - the "convenient interface for using the LP solver" could be factored
#   out into a separate object.

import argparse
import math
import pdb
import sys

import itertools
import more_itertools
import numpy as np
import scipy.optimize
import scipy.sparse
# note that comb() returns a float by default; for loop bounds, it needs
# the "exact=True" option, so that it returns an int
from scipy.special import comb
from scipy.stats import hypergeom

def cliques_left_after_zeroing(clique_set, edge):
    """Finds cliques which are left after zeroing out an edge.

    clique_set: a set of cliques
    edge: an edge (as a two-element set)
    Returns: the cliques in clique_set which aren't hit by that edge.
    """
    return frozenset([c for c in clique_set if not edge < c])

class LatticeRankBound:
    """Computes a bound on the rank of functions for finding cliques.

    Here, we try to emphasize code simplicity over speed.
    This includes the full lattice! So it has a ton of variables.
    (Hopefully, these can eventually be grouped somehow.)
    """
    def __init__(self, n, k):
        """ Constructor.

        n: number of vertices
        k: size of clique we're looking for
        """
        # problem size
        self.n = n
        self.k = k
        # list of all the possible cliques
        self.all_cliques = [frozenset(c)
            for c in itertools.combinations(range(n), k)]
        # the variables are alll possible sets of cliques: this maps them
        # to a variable index in the LP
        self.var_index = {}
        for s in more_itertools.powerset(self.all_cliques):
            self.var_index[frozenset(s)] = len(self.var_index)

        # These store the constraints:
        # A: a list of lists of (A,i,j) entries (which go into a sparse matrix)
        # b: a list of numbers
        # the inequalities (note that the LP solver expects upper bounds)
        self.A_ub = []
        self.b_ub = []
        # the equalities, stored similarly
        self.A_eq = []
        self.b_eq = []

    def add_constraint(self, A, op, b):
        """Adds one row to the constraints.

        A: a list of (index, coefficient) pairs, where "index" is
            a key (of any hashable Python type) in var_index
        op: either '<', '=', or '>': this is the type of constraint
            ??? can we treat '>' the same as '>='?
        b: the corresponding bound (a float)
        Side effects: adds the constraint
        """
        # print(str(A) + ' ' + op + ' ' + str(b))
        # converts from "list of coefficient" to a row of A
        def get_coefs(i, negate): 
            # this is arguably pushing limits for a list comprehension...
            return [(-a if negate else a, i, self.var_index[k])
                for (k,a) in A]
        # add whichever kind of constraint
        if op == '<':
            i = len(self.b_ub)
            self.A_ub += get_coefs(i, False)
            self.b_ub += [b]
            return
        if op == '=':
            i = len(self.b_eq)
            self.A_eq += get_coefs(i, False)
            self.b_eq += [b]
            return
        if op == '>':
            i = len(self.b_ub)
            self.A_ub += get_coefs(i, True)
            self.b_ub += [-b]
            return

    def add_average_of_all_sets_constraint(self):
        """Adds constraint on the average of all the ranks."""
        num_sets = len(self.var_index.keys())
        # this is the average, over all of the sets
        A = [(s, 1.0/num_sets) for s in self.var_index.keys()]
        b = (num_sets-1) / 2.0
        self.add_constraint(A, '>', b)

    def add_edge_zeroing_constraints(self):
        """Add constraints that "zeroing an edge removes some gates".

        This is arguably a relatively simple constraint: it's
        basically encoding the edges in Z(X,Y).
        """
        # loop through the edges
        for edge in itertools.combination(range(self.n), 2):
            # loop through sets of cliques
            for Y in self.var_index.keys():
                # find effect of zeroing out that edge
                X = cliques_left_after_zeroing(Y, frozenset(edge))
                # is this smaller? (that is, did the edge "hit" B?
                if X < Y:
                    # constrain the "sets above" X to have more gates
                    # (and so has higher rank)
                    self.add_constraint([(Y, 1.0), (X, -1.0)], '>', 0)

    def add_higher_set_constraints(self):
        """Constrains sets containing an edge e to be above sets without e.

        Let X be a set of cliques, none of which include an edge e.
        Then if Y is X + (any set of cliques including e), |C(Y)| >= |C(X)|.
        """
        # FIXME not yet implemented
        # loop through the edges
        for e in itertools.combination(range(self.n), 2):
            # get the sets of cliques which include e
            B_sets = [s for s in self.var_index.keys() if contains_edge(s, e)]
            # loop through sets of cliques
            for A in self.var_index.keys():
                # if this set includes e, skip it
                if cliques_left_after_zeroing(A, e) < A:
                    continue
                # constrain the "sets above" A to be higher
                pass   # FIXME

    def solve(self):
        """Solves the linear system.
        
        Note that by default, the solver constrains all x >= 0,
        so we don't add that constraint.
        Returns: the LP result
        """
        # utility to convert entries to a sparse array
        def sparse_array_from_entries(A):
            # gets i'th element of all the tuples
            def ith(i):
                return [a[i] for a in A]
            return scipy.sparse.coo_array( (ith(0), (ith(1), ith(2))) )
        # convert A and b to np.array objects
        A_ub = sparse_array_from_entries(self.A_ub)
        b_ub = np.array(self.b_ub)
        A_eq = sparse_array_from_entries(self.A_eq)
        b_eq = np.array(self.b_eq)

        # the objective function: how low can the rank of finding
        # all the cliques be?
        c = np.zeros(len(self.var_index))
        c[ self.var_index[frozenset(self.all_cliques)] ] = 1
        # solve
        # ??? Is there a way to tell the solver that this is sparse?
        # (It's detecting this, but that throws a warning.)
        r = scipy.optimize.linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq)
        # FIXME deal with this failing
       
        # FIXME return average for each level (not just the bound for CLIQUE)
        pdb.set_trace()

        return x

    def get_bound(self, include_upper_bound):
        """Gets the bound.

        include_upper_bound: if True, include the upper bound
        Returns: a 2-D NumPy array, with axes "# cliques zeroed"
            and "# cliques remaining"
        """
        self.add_total_cliques_equality_constraints()
        self.add_total_cliques_counting_bound_constraints()
        self.add_edge_zeroing_constraints()
        # possibly include the upper bound
        if include_upper_bound:
            self.add_upper_bound_constraints()
        x = self.solve()
        return x

if __name__ == '__main__':
    # gate_bound_smoke_test()

    parser = argparse.ArgumentParser(
        description='Attempt at bound on finding some number of cliques.')
    parser.add_argument('n', type=int,
        help='number of vertices in input graph')
    parser.add_argument('k', type=int,
        help='number of vertices in cliques to find')
    parser.add_argument('--include-upper-bound', action='store_true',
        help='include the upper bound constraint')
    parser.add_argument('--plot-surface', action='store_true',
        help='plot the bounds as a surface, with and without the upper bound')
    args = parser.parse_args()

    # possibly plot the surfaces
    if args.plot_surface:
        # plot the surfaces
        plot_bound_surfaces(args.n, args.k, 'bound_surfaces_')
        sys.exit(0)

    # otherwise, print the bound
    lp = LpBound(args.n, args.k)
    x = lp.get_bound(args.include_upper_bound)
    bound = x[ lp.max_cliques_zeroed, lp.max_cliques_remaining ]
    print(np.round(bound, 4))

