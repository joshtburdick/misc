#!/usr/bin/env python3
# Convenient interface (hopefully) for use with PuLP.

import fractions
import math
import os
import pdb
import re
import sys
import subprocess
import tempfile

import functools
import numpy as np

import pulp

class PuLP_Helper:
    """Wrapper class for PuLP solver, providing convenient variable names.
    """
    def __init__(self, var_names):
        """ Constructor.
    
        var_names: names of the variables
        """
        # mapping from variable name (which needn't be a string)
        # to PuLP variable object
        self.vars = {}
        # the variables are assumed to all be >= 0
        for v in var_names:
            self.vars[v] = pulp.LpVariable(
                self.get_parseable_name(v), 0)
        # the problem
        self.prob = pulp.LpProblem("sidsproblem", pulp.LpMinimize)

    def get_parseable_name(self, var_name):
        """Makes a 'parseable' version of a variable name."""
        s = str(var_name)
        s = re.sub("[ ']", "", s)
        s = re.sub("[\\.,\\(\\)]", "_", s)
        return "x_" + s

    def add_constraint(self, A, op, b):
        """Adds one row to the problem.

        A: a list of (x, a) pairs, where:
            a is the coefficient
            x is the variable, which is a key in self.vars
        (XXX arguably "(a, x)" would make sense, but I used "(x, a)"
            elsewhere in the code, so using that)
        op: either '<=', '==', or '>=': the type of constraint
            (??? change "=" to "==" ?)
        b: the corresponding bound
        Side effects: adds the constraint
        """
        if op not in ["<=", "=", ">="]:
            raise ValueError(f"unknown operator: {op}")
        # convert coefficients to format PuLP expects
        A_as_expr = pulp.lpSum([a * self.vars[x] for (x, a) in A])
        if op == "<=":
            self.prob += A_as_expr <= b
        if op == "=":
            self.prob += A_as_expr == b
        if op == ">=":
            self.prob += A_as_expr >= b

    def solve_1(self, var_to_minimize):
        """Solves the linear system, for one variable.

        This assumes all variables are >= 0.
        FIXME add option to make all variables integers?
        """
        self.prob += self.vars[var_to_minimize]
        # for debugging
        # self.prob.writeLP("./bound.lp")
        r = self.prob.solve(pulp.GLPK(options=['--exact']))
        # problem had a solution
        if r == 1:
            return self.vars[var_to_minimize].varValue
        # problem was infeasible, or something else went wrong
        else:
            return None

    def solve(self, var_to_minimize, bounds=None):
        """Solves the linear system.

        objective: what to minimize (as a Numpy vector)
        var_to_minimize: name of the variable to minimize
        bounds: bounds on individual variables
        Returns: a dict, indexed by variable name, of
            all the variables, at the lower bound
        """
        self.prob += self.vars[var_to_minimize]
        self.prob.writeLP("./bound.lp")
        r = self.prob.solve(pulp.GLPK(options=['--exact']))
        # r = self.prob.solve(pulp.GLPK())
        # problem had a solution
        print(f"Result r = {r}")
        if r == 1:
            opt = {x: self.vars[x].varValue
                for x in self.vars}
            return opt
        # problem was infeasible, or something else went wrong
        else:
            return None

