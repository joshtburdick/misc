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
    def __init__(self, var_names, verbosity=1):
        """ Constructor.
   
        var_names: names of the variables
        """
        self.verbosity = verbosity
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
        i = len(self.var_index)
        self.var_index[f"__slack{self.n_slack_vars}"] = i
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

        # following Burke textbook, for equalities, add "<=" and ">=" ?
        # if op == "=":
        #     self.add_constraint(A, "<=", b)
        #     self.add_constraint(A, ">=", b)
        #     return

        # at this point, op is either "<=" or ">="
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


    def solve(self, var_to_optimize, bounds=None, minimize=True):
        """Solves the linear system, optimizing for one variable.

        var_to_optimize: name of var to minimize or maximize
        bounds: the bounds (currently ignored)
        minimize: if True, minimize; otherwise maximize
        """
        # the (sparse) objective function
        c = [(var_to_optimize, (-1 if minimize else 1))]
        return self.solve_1(c)


    def solve_1(self, objective):
        """Solves the linear system.

        objective: what to minimize (as a Numpy vector)
        var_to_minimize: name of the variable to minimize
        Returns: the vector of the minizing solution (if found),
            as a Numpy vector.  FIXME: return a dict, indexed by variable name, of
            all the variables, at the lower bound?
        """
        # convert from var. names to indices
        c = {self.var_index[name]: fractions.Fraction(x)
            for (name, x) in objective }
        if self.verbosity >= 1:
            print(f"var_index = {self.var_index}")
            print(f"A =")
            for r in self.A.items():
                print(r)
            print(f"b =\n{self.b}")
            print(f"c =\n{c}")
        # run the simplex algorithm
        # t, s, v = exactsimplex.sparse.simplex(c, self.A, self.b)
        t, s, v = exactsimplex.sparse.simplex_two_phase(c, self.A, self.b)
        # pdb.set_trace()
        # FIXME should check for errors
        opt_vec = {var: (s[i] if i in s else 0)
            for (var, i) in self.var_index.items()}
        return opt_vec


