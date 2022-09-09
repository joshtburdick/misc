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

    def add_average_of_levels_constraint(self, max_cliques, constraint_type):
        """Adds constraint on the average of some levels.

        max_cliques: this will average sets including from 0 to
            this number (inclusive) of cliques
        constraint_type: either '>' or '='.
            Note that '=' is only valid if
            max_cliques == the total number of cliques.
        Side effects: adds a constraint on the average of levels from
            0 to max_cliques (inclusive).
        """
        # get all the sets which are small enough
        clique_sets = [s for s in self.var_index.keys()
            if len(s) <= max_cliques]
        num_sets = len(clique_sets)
        # this is the average, over all of the sets
        A = [(s, 1.0/num_sets) for s in clique_sets]
        # the average of their ranks (using 0-based numbering)
        b = (num_sets-1) / 2.0
        # add the constraint
        self.add_constraint(A, constraint_type, b)

    def add_edge_zeroing_constraints(self):
        """Add constraints that "zeroing an edge removes some gates".

        This is arguably a relatively simple constraint: it's
        basically encoding the edges in Z(X,Y).
        """
        # loop through the edges
        for edge in itertools.combinations(range(self.n), 2):
            # loop through sets of cliques
            for Y in self.var_index.keys():
                # find effect of zeroing out that edge
                X = cliques_left_after_zeroing(Y, frozenset(edge))
                # is this smaller? (that is, did the edge "hit" B?
                if X < Y:
                    # constrain the "sets above" X to have more gates
                    # (and so has higher rank)
                    self.add_constraint([(Y, 1.0), (X, -1.0)], '>', 1)

    def add_higher_set_constraints(self):
        """Constrains sets containing an edge e to be above sets without e.

        Let A be a set of cliques, none of which include an edge e.
        Let B be any set of cliques including e.
        Then E[rank(A+B)] = E[rank(A)] + (|B|+1)/2 .

        This doesn't appear to help much. Also, it's more complicated.
        """
        # loop through the edges
        for e in itertools.combinations(range(self.n), 2):
            # get the sets of cliques which include e
            B_cliques = [c for c in self.all_cliques if c > frozenset(e)]
            # include all nonempty combinations of these
            B_sets = [frozenset(s) for s in more_itertools.powerset(B_cliques)
                if len(s) > 0]
            # loop through sets of cliques
            for A in self.var_index.keys():
                # if this set includes e, skip it
                if cliques_left_after_zeroing(A, frozenset(e)) < A:
                    continue
                # constrain the average of A and "sets above" A
                # to be higher than A
                B_set_constraints = [(A.union(B), 1./len(B_sets))
                    for B in B_sets]
                self.add_constraint(B_set_constraints + [(A, -1.)],
                    '>',
                    # since these are all distinct, and "higher" than A, this
                    # is the average of numbers from 1 to |B_sets|, inclusive
                    (len(B_sets)+1.) / 2)

    def solve(self):
        """Solves the linear system.

        Returns: a vector, indexed by the number of cliques, of the
            lower bound on the average for that level
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
        # possibly add equality constraints
        A_eq = None
        b_eq = None
        if self.A_eq:
            A_eq = sparse_array_from_entries(self.A_eq)
            b_eq = np.array(self.b_eq)
        # the objective function: how low can the rank of finding
        # all the cliques be?
        c = np.zeros(len(self.var_index))
        c[ self.var_index[frozenset(self.all_cliques)] ] = 1
        # This is the maximum rank: if we are ranking all of the functions,
        # it can't go higher than the number of subsets of cliques.
        # (The "-1" is because we're numbering ranks starting at 0.)
        max_rank = len(self.var_index) - 1
        # solve
        # ??? Is there a way to tell the solver that this is sparse?
        # (It's detecting this, but that throws a warning.)
        r = scipy.optimize.linprog(c, A_ub=A_ub, b_ub=b_ub,
            A_eq=A_eq, b_eq=b_eq,
            # we include these bounds, although they don't seem to help
            # all that much
            bounds = (0, max_rank))
        # FIXME deal with this failing
        # FIXME move everything after this to get_bound() ?
        # get average for each "level" (not just the bound for CLIQUE)
        total_rank = np.zeros(len(self.all_cliques)+1)
        num_functions = np.zeros(len(self.all_cliques)+1)
        # loop through the cliques, and average the ranks
        for (s, i) in self.var_index.items():
            num_cliques = len(s)
            total_rank[num_cliques] += r.x[i]
            num_functions[num_cliques] += 1
        rank_bound = total_rank / num_functions
        # pdb.set_trace()
        return rank_bound

    def get_bound(self, include_expectation_lower_bound,
            include_zeroing_bound, include_expectation_equality_bound):
        """Gets the bound.

        include_expectation_lower_bound: if True, include the lower bound
            on expected values of levels
        include_zeroing_bound: if True, include bound from zeroing out
            each edge
        include_expectation_equality_bound: if True, include the exact
            bound on expected value of _all_ the levels
        Returns: a 1-D Numpy array, whose i'th element is a
            lower bound on the rank of finding i cliques
        """
        # include bounds on expected value of all the levels?
        if include_expectation_lower_bound:
            # loop through the levels of numbers of cliques
            # n.b.: adding all these levels seems to help more than
            # only adding the >= constraint for "all the cliques"
            for i in range(len(self.all_cliques)+1): 
                 self.add_average_of_levels_constraint(i, '>')
        # include bounds from zeroing edges?
        if include_zeroing_bound:
            self.add_edge_zeroing_constraints()
        # include exact bound on average of _all_ the levels?
        if include_expectation_equality_bound:
            self.add_average_of_levels_constraint(len(self.all_cliques)+1, '=')

        # N.B.: this doesn't seem to help much. But tossing it in anyway...
        self.add_higher_set_constraints()
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
    args = parser.parse_args()
    # get the bound
    lp = LatticeRankBound(args.n, args.k)
    x = lp.get_bound(True, True, True)
    # pdb.set_trace()
    print(x)
    # for now, just printing the bound for CLIQUE
    bound = x[ len(lp.all_cliques) ]
    print(np.round(bound, 4))

