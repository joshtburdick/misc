#!/usr/bin/env python3
# Convenient interface (hopefully) for using an SCIPopt.
# ??? is this actually arbitrary precision?
# ??? see also https://pypi.org/project/flexible-lp/
# FIXME
# - support other LP solvers?
#   - if that happens, rename this?

import collections
import fractions
import io
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

    def parse_scip_output(self, solution_file):
        """Parses SCIP's output."""
        with open(solution_file, "rt") as f:
            line = f.readline()
            if not re.match(r"solution status: optimal solution found", line):
                return None
            line = f.readline()
            if not re.match(r"objective value:\s+(\w+)", line):
                return None
            # remaining lines are mostly of the form:
            # x__E_10_          1   (obj:1)
            # this dict will be indexed by "parseable variable name"
            x_by_parseable_name = collections.defaultdict(float)
            for line in f: 
                m = re.match(r"(\w+)\s+(\S+) \t", line)
                x_by_parseable_name[m.group(1)] = float(m.group(2))
            x_by_name = {name: x_by_parseable_name[parseable_name]
                for name, parseable_name in self.var_name.items()
                if parseable_name in x_by_parseable_name
            }
            return x_by_name

    def solve_1(self, objective, bounds=None):
        """Solves the linear system, for one variable.

        objective: what to minimize (as a Numpy vector)
        var_to_minimize: name of the variable to minimize
        bounds: bounds on individual variables
        Returns: the optimal value of the objective function,
            or None if the problem was infeasible.
        """
        # may not use this...
        raise NotImplementedError()

    def solve(self, var_to_minimize, lp_filename="./bound.lp"):
        """Solves the linear system.

        This assumes all variables are >= 0.
        FIXME add option to make all variables integers?
        """
        # construct the linear program
        newline = "\n"     # ??? possibly python3.8 workaround?
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
            script = f"""
read {lp_filename}
optimize
write solution bound_opt.txt
quit
"""
            scip_output = subprocess.check_output(["scip"], text=True, input=script)
            return self.parse_scip_output("bound_opt.txt") 

            ##### rest of this goes away...
            # parse line of the syntax, e.g.:
            # Primal Bound     : +1.22500000000000e+02   (in run 1, after 1 nodes, 0.02 seconds, depth 2, found by <locks>)
#            m = re.search(".*Primal Bound\s+:[ ]*([^ ]+)[ ]+\(in run", scip_output)
#  First Solution   : +0.00000000000000e+00   (in run 1, after 0 nodes, 0.01 seconds, depth 0, found by <relaxation>)
            m = re.search(r".*First Solution\s*:[ ]*([^ ]+)[ ]*\(in run.*", scip_output)
            if not m:
                return None
            bound = float(m.group(1))
        return bound


