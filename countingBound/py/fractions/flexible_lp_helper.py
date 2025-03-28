#!/usr/bin/env python3
# Convenient interface (hopefully) for using the `flexible-lp` library.

import math
import pdb
import sys

import fractions
import numpy as np
import scipy.optimize
import scipy.sparse

import flexible_lp.simplex

class Flexible_LP_Helper:
    """Wrapper class for LP solver, providing convenient variable names.

    Uses the flexible-lp library.
    """
    def __init__(self, var_names):
        """ Constructor.
    
        var_names: names of the variables
        """
        # mapping from variable name to (numerical) index
        self.var_index = {}
        for v in var_names:
            self.var_index[v] = len(self.var_index)
        # These store the constraints, each as a dict with keys:
        # A: a list of lists of (A,i,j) entries (which go into a sparse matrix)
        # op: either ">=", "=", or ">="
        # b: a list of numbers
        self.constraints = []

    def add_constraint(self, A, op, b):
        """Adds one row to the constraints.

        A: a list of (index, coefficient) pairs, where "index" is
            a key (of any hashable Python type) in var_index
        op: either '<', '<=', '=', '>=', or '>': the type of constraint
        b: the corresponding bound
        Side effects: adds the constraint
        """
        # print(str(A) + ' ' + op + ' ' + str(b))
        def negate_coefs(A):
            """Negates coefficients in a list."""
            return [(a,-x) for (a,x) in A]
        # for some operators, simply add the constraint
        if op in ["<=", "=", ">="]:
            self.constraints.append({
                "A": A,
                "op": op,
                "b": b
            })
            return
        # for other operators, need negation (although I'm not sure
        # these will be used)
        if op in [">", "<"]:
            self.constraints.append({
                "A": negate_coefficients(A),
                "op": ">=" if op=="<" else "<=",
                "b": -b
            })
            return
        raise ValueError(f"Unknown operator: {op}")

    def solve(self, var_to_minimize, bounds=None):
        """Solves the linear system, minimizing one variable.

        Note that this ignores the bounds.
        """
        # construct a vector with a 1. only at the variable
        # to minimize (and 0. everywhere else)
        c = [fractions.Fraction(0)] * len(self.var_index)
        c[ self.var_index[var_to_minimize] ] = fractions.Fraction(1)
        opt_val, opt_vec = self.solve_1(c)
        # check whether problem was infeasible
        if opt_val is None:
            return None
        # FIXME: convert the result to a dict, and return opt_vec?
        opt_vec = {var: opt_vec[i]
            for (var, i) in self.var_index.items()}
        return opt_val

    def solve_1(self, objective):
        """Solves the linear system.

        objective: what to minimize (as a Numpy vector)
        var_to_minimize: name of the variable to minimize
        Returns: the vector of the minizing solution (if found),
            as a Numpy vector.
        FIXME: return a dict, indexed by variable name, of
            all the variables, at the lower bound?
        """
        def coefs_to_list(A):
            """Converts a sparse list of coefficients to a list of numbers."""
            A_dict = {self.var_index[x]: a for (x, a) in A}
            A_array = [fractions.Fraction(A_dict[i]) if i in A_dict else fractions.Fraction(0)
                for i in range(len(self.var_index))]
            return A_array
        def get_constraints_with_op(op):
            """Gets constraints with a given operator."""
            constraints = [c for c in self.constraints
                if c["op"] == op]
            return ([coefs_to_list(c["A"]) for c in constraints],
                [fractions.Fraction(c["b"]) for c in constraints])
        # gather the constraints
        (A_g, b_g) = get_constraints_with_op(">=")
        (A_e, b_e) = get_constraints_with_op("=")
        (A_l, b_l) = get_constraints_with_op("<=")
        # pdb.set_trace()
        # call the solver
        (opt_val, opt_vec) = flexible_lp.simplex.linprog(objective,
            A_g=A_g, b_g=b_g, A_e=A_e, b_e=b_e, A_l=A_l, b_l=b_l,
            value_map=fractions.Fraction)
        return (opt_val, opt_vec)

