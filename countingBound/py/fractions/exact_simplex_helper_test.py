# Some small tests.

import unittest

import exact_simplex_helper

class TestExactSimplexHelper(unittest.TestCase):

    def test1(self):
        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2", "y"])
        h.add_constraint([("x1", 1), ("x2", 5)], ">=", -10)
        h.add_constraint([("x1", 5), ("x2", 1)], ">=", -10)
        h.add_constraint([("x1", 1), ("x2", 1), ("y", -1)], "=", 0)
        r = h.solve("y")
        print(r)




