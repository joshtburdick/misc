#!/usr/bin/env python3
# Convenient interface (hopefully) for using an SCIPopt.
# ??? is this actually arbitrary precision?
# ??? see also https://pypi.org/project/flexible-lp/
# FIXME
# - support other LP solvers?
#   - if that happens, rename this?

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
import scipy.optimize
import scipy.sparse

class SCIP_Helper:
    """Wrapper class for SCIP solver, providing convenient variable names.

    Also deals with fractional coefficients.
    """
    def __init__(self, var_names):
        """ Constructor.
    
        var_names: names of the variables
        """
        # mapping from variable name (which needn't be a string)
        # to "parseable name" (which can be parsed by SCIPopt)
        self.var_name = {}
        for v in var_names:
            self.var_name[v] = self.get_parseable_name(v)
        # this stores the constraints, as strings to pass to SCIPopt
        self.constraints = []

    def get_parseable_name(self, var_name):
        """Makes a 'parseable' version of a variable name."""
        s = str(var_name)
        s = re.sub("[ ']", "", s)
        s = re.sub("[\\.,\\(\\)]", "_", s)
        return "x_" + s

    def as_int(self, x):
        """Utility converting a fraction to an integer."""
        x = fractions.Fraction(x)
        if x.is_integer() and x.denominator==1:
            return x.numerator
        return None

    def add_constraint(self, A, op, b):
        """Adds one row to the constraints.

        In this case, this involves:
        - scaling the constraint to get rid of fractions, and
        - light formatting into LP format.
        A: a list of (index, coefficient) pairs, where "index" is
            a key (of any hashable Python type) in var_index
        op: either '<', '=', or '>': this is the type of constraint
            ??? can we treat '>' the same as '>='?
        b: the corresponding bound (a float)
        Side effects: adds the constraint
        """
        if op not in ["<", "<=", "=", ">=", ">"]:
            raise ValueError(f"unknown operator: {op}")
        # get least common multiple of denominators
        coefs = [fractions.Fraction(a) for (_,a) in A]
        denominators = [fractions.Fraction(x).denominator 
            for x in coefs + [b]]
        m = functools.reduce(math.lcm, denominators)
        # get constraint, as strings, with fractions cancelled out by
        # multiplying by LCM
        A_as_strings = [str(self.as_int(m*a)) + " " + self.var_name[name]
                        for (name, a) in A]
        constraint = " + ".join(A_as_strings) + " " + op + " " + str(self.as_int(m*b))
        # print(constraint)
        self.constraints.append(constraint)

    def solve(self, var_to_minimize, lp_filename="./bound.lp"):
        """Solves the linear system, for one variable.

        This assumes all variables are >= 0.
        FIXME add option to make all variables integers?

        """
        # construct the linear program
        newline = "\n"   # ??? possibly python3.8 workaround?
        lp_string = f"""
Minimize
 obj: {self.var_name[var_to_minimize]}
Subject To
{ newline.join(self.constraints) }
Bounds
{ newline.join([f'0 <= {x}' for x in self.var_name.values()])}
End
"""
        # XXX for debugging
        # with tempfile.NamedTemporaryFile("wt", suffix=".lp") as lp_file:
        os.makedirs(os.path.dirname(lp_filename), exist_ok=True)
        with open(lp_filename, "wt") as lp_file:
            lp_file.write(lp_string)
            lp_file.flush()
            scip_output = subprocess.check_output([
                "scip", "-f", lp_file.name], text=True)
            # pdb.set_trace()
            # print(scip_output)
            # pdb.set_trace()
            # parse line of the syntax, e.g.:
            # Primal Bound     : +1.22500000000000e+02   (in run 1, after 1 nodes, 0.02 seconds, depth 2, found by <locks>)
#            m = re.search(".*Primal Bound\s+:[ ]*([^ ]+)[ ]+\(in run", scip_output)
#  First Solution   : +0.00000000000000e+00   (in run 1, after 0 nodes, 0.01 seconds, depth 0, found by <relaxation>)
            m = re.search(".*First Solution\s*:[ ]*([^ ]+)[ ]*\(in run.*", scip_output)
            if not m:
                return None
            print(m.group(1))
            bound = float(m.group(1))
        return bound

    def solve_1(self, objective, bounds=None):
        """Solves the linear system.

        objective: what to minimize (as a Numpy vector)
        var_to_minimize: name of the variable to minimize
        bounds: bounds on individual variables
        Returns: the vector of the minizing solution (if found),
            as a Numpy vector.
        FIXME: return a dict, indexed by variable name, of
            all the variables, at the lower bound?
        """
        # for now, not bothering with this
        raise NotImplementedError()

