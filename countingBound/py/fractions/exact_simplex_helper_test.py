# Some small tests.

import unittest

import exact_simplex_helper

class TestExactSimplexHelper(unittest.TestCase):

    def test1(self):
        """Example from explication."""
        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2", "y"])
        # since this only supports minimizing a single variable,
        # including it with an equality constraint
        h.add_constraint([("x1", 3), ("x2", 2), ("y", -1)], "=", 0)
        h.add_constraint([("x1", 1), ("x2", 2)], "<=", 4)
        h.add_constraint([("x1", 1), ("x2", -1)], "<=", 1)
        r = h.solve("y")
        print(r)


    def test2(self):
        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2", "y"])
        h.add_constraint([("x1", 1), ("x2", 5)], ">=", 10)
        h.add_constraint([("x1", 5), ("x2", 1)], ">=", 10)
        h.add_constraint([("x1", 1), ("x2", 1), ("y", -1)], "=", 0)
        r = h.solve("y")
        print(r)




