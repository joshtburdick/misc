#!/usr/bin/env python3
# Convenient interface (hopefully) for using the `simplex_algorithm` library
# (from Jules).

import math
import pdb
import sys

import fractions
import numpy as np
import scipy.optimize
import scipy.sparse

from simplex_algorithm.simplex import SimplexSolver

class SimplexAlgorithmHelper:
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
        # A stores the rows of the tableau, as a sparse 2-D array,
        # represented as a list of rows; each row is represented by a
        # a dict of fractions.Fraction objects .
        self.A = []
        # b contains the last column, as a list
        self.b = []

    def add_constraint(self, A, op, b):
        """Adds one row to the constraints.

        A: a list of (index, coefficient) pairs, where "index" is
            a key (of any hashable Python type) in var_index
        op: either '<=', '=', or '>=': the type of constraint
        b: the corresponding bound
        Side effects: adds the constraint
        """
        b = fractions.Fraction(b)
        if op not in ["<=", "=", ">="]:
            raise ValueError(f"Unknown operator: {op}")
        # print(str(A) + ' ' + op + ' ' + str(b))

        # For equalities, add both "<=" and ">="
        if op == "=":
            self.add_constraint(A, "<=", b)
            self.add_constraint(A, ">=", b)
            return

        # At this point, op is either "<=" or ">=";
        # if it's ">=", then negate A's entries, and b.
        sign = -1 if op == ">=" else 1

        self.A.append({self.var_index[name]: sign * fractions.Fraction(x)
            for (name, x) in A})
        self.b.append( sign * fractions.Fraction(b) )


    def dict_to_list(self, x):
        """Converts a dict (keyed by integers) to a list.

        Missing values are converted to 0's.
        This is used both for the rows of A, and the objective function.
        (Ideally, we wouldn't be constructing a list with lots of 0's.
        If the solver could use a sparse input, this wouldn't be needed.)
        """
        def lookup(i):
            try:
                return x[i]
            except KeyError:
                return fractions.Fraction(0)
        return [lookup(i) for i in range(len(self.varnames))]

    def solve(self, var_to_optimize, minimize=True):
        """Solves the linear system, optimizing for one variable.

        var_to_optimize: name of var to minimize or maximize
        minimize: if True, minimize; otherwise maximize
        """
        # the (sparse) objective function
        c = [(var_to_optimize, (-1 if minimize else 1))]
        return self.solve_1(c)

    def solve_1(self, objective):
        """Solves the linear system.

        objective: what to minimize (as a Numpy vector)
        var_to_minimize: name of the variable to minimize
        Returns: a dict, indexed by variable name, of
            all the variables, at the lower bound
        """
        # convert from var. names to indices
        c = {self.var_index[name]: fractions.Fraction(x)
            for (name, x) in objective }
        if self.verbosity >= 1:
            print(f"var_index = {self.var_index}")
            print(f"A =")
            for r in self.A:
                print(r)
            print(f"b =\n{self.b}")
            print(f"c =\n{c}")
        # run the simplex algorithm
        solver = SimplexSolver(c, self.A, self.b, sparse_input=True)
        status = solver.solve(verbose=False)
        if status == "optimal":
            solution = solver.get_solution()
            # FIXME: Currently, the solver returns a dict, keyed by
            # "x1", "x2", "x3", etc. If it returned a dict, indexed by
            # 0-based variable index, this would be simpler.
            opt_vec = {var: (float(solution[f"x{i+1}"]) if f"x{i+1}" in solution else 0)
                for (var, i) in self.var_index.items()}
            # pdb.set_trace()
            return opt_vec
        # FIXME should check for errors
        else:
            print(f"solver returned status = {status}")
            return None

