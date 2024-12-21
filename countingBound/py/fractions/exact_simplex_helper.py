#!/usr/bin/env python3
# Convenient interface (hopefully) for using the `exactsimplex` library.

import math
import pdb
import sys

import fractions
import numpy as np
import scipy.optimize
import scipy.sparse

import exactsimplex.sparse

class ExactSimplexHelper:
    """Wrapper class for LP solver, providing convenient variable names.

    This also handles converting the LP to canonical form.
    """
    def __init__(self, var_names):
        """ Constructor.
   
        var_names: names of the variables
        """
        # mapping from variable name to (numerical) index
        self.var_index = {}
        for v in var_names:
            self.var_index[v] = len(self.var_index)
        # A stores the rows of the tableau, as a sparse 2-D array.
        # The outer dict contains the rows; each row is represented by a
        # a dict of fractions.Fraction objects .
        self.A = {}
        # b, similarly, contains the last column
        self.b = {}
        # The number of slack variables, which will be prefixed with __slack
        self.n_slack_vars = 0

    def allocate_slack_var(self):
        """Allocates a new slack variable."""
        i = self.n_slack_vars
        self.var_index[f"__slack{i}"] = len(self.var_index)
        self.n_slack_vars += 1
        return i

    def add_constraint(self, A, op, b):
        """Adds one row to the constraints.

        A: a list of (index, coefficient) pairs, where "index" is
            a key (of any hashable Python type) in var_index
        op: either '<=', '=', or '>=': the type of constraint
        b: the corresponding bound
        Side effects: adds the constraint
        """
        if op not in ["<=", "=", ">="]:
            raise ValueError(f"Unknown operator: {op}")
        # print(str(A) + ' ' + op + ' ' + str(b))
        # build up the row, starting with A 
        row = {self.var_index[name]:
            fractions.Fraction(x) for (name, x) in A}
        # add slack vars if necessary
        if op == "<=":
            row[ self.allocate_slack_var() ] = 1
        if op == ">=":
            row[ self.allocate_slack_var() ] = -1
        # add row to A, and add to b
        i = len(self.A)
        self.A[i] = row
        self.b[i] = fractions.Fraction(b)

    def solve(self, var_to_minimize, bounds=None):
        """Solves the linear system, minimizing one variable.

        Note that this ignores the bounds.
        """
        # the (sparse) objective function; since we're
        # minimizing, put -1 for the variable in question
        c = { self.var_index[var_to_minimize]: -1 }
        pdb.set_trace()
        # run the simplex algorithm
        t, s, v = exactsimplex.sparse.simplex(c, self.A, self.b)
        # FIXME simplex should check for infeasible problems
        opt_vec = {var: (s[i] if i in s else 0)
            for (var, i) in self.var_index.items()}
        return opt_vec

    def solve_1(self, objective):
        """Solves the linear system.

        objective: what to minimize (as a Numpy vector)
        var_to_minimize: name of the variable to minimize
        Returns: the vector of the minizing solution (if found),
            as a Numpy vector.
        FIXME: return a dict, indexed by variable name, of
            all the variables, at the lower bound?
        """
        raise NotImplementedError()

