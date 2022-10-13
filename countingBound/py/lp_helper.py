#!/usr/bin/env python3
# Convenient interface (hopefully) for using the LP solver.

import math
import pdb

import numpy as np
import scipy.optimize
import scipy.sparse

class LP_Helper:
    """Wrapper class for LP solver, providing convenient names.

    """
    def __init__(self, var_names):
        """ Constructor.
    
        var_names: names of the variables
        """
        # mapping from variable name to index
        self.var_index = {}
        for v in var_names:
            self.var_index[v] = len(self.var_index)
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
            # express "Ax > b" as "-Ax < -b"
            self.A_ub += get_coefs(i, True)
            self.b_ub += [-b]
            return

    def solve(self, var_to_minimize, bounds=None):
        """Solves the linear system.

        var_to_minimize: name of the variable to minimize
        bounds: bounds on individual variables
        Returns: a dict, indexed by variable name, of
            the lower bound.
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
        c = np.zeros(len(self.var_index))
        c[ self.var_index[var_to_minimize] ] = 1
        # solve
        # ??? Is there a way to tell the solver that this is sparse?
        # (It's detecting this, but that throws a warning.)
        r = scipy.optimize.linprog(c, A_ub=A_ub, b_ub=b_ub,
            A_eq=A_eq, b_eq=b_eq,
            # we include these bounds, although they don't seem to help
            # all that much
            bounds = bounds)
        # pdb.set_trace()
        # FIXME deal with this failing
        # convert the bound to a dict
        bound = {var: r.x[i]
            for (var, i) in self.var_index.items()}
        return bound
