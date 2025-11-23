# Tests, comparing with scipy's implementation.


import fractions

import scipy.optimize
import unittest

import exact_simplex_helper

class TestExactSimplexHelper2(unittest.TestCase):


    def test1(self):
        """Example from explication."""

        FIXME

        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2"])
        h.add_constraint([("x1", 1), ("x2", 2)], "<=", 4)
        h.add_constraint([("x1", 1), ("x2", -1)], "<=", 1)
        r = h.solve_1([("x1", 3), ("x2", 2)])
        self.assertEqual(r["x1"], 2)
        self.assertEqual(r["x2"], 1)
        print(r)





