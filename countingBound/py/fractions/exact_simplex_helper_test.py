# Some small tests.

import unittest

import exact_simplex_helper

class TestExactSimplexHelper(unittest.TestCase):

    def test1(self):
        """Example from explication."""
        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2"])
        h.add_constraint([("x1", 1), ("x2", 2)], "<=", 4)
        h.add_constraint([("x1", 1), ("x2", -1)], "<=", 1)
        r = h.solve_1([("x1", 3), ("x2", 2)])
        print(r)


    def skip_test2(self):
        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2", "y"])
        h.add_constraint([("x1", 1), ("x2", 5)], "<=", 10)
        h.add_constraint([("x1", 5), ("x2", 1)], "<=", 10)
        h.add_constraint([("x1", 1), ("x2", 1), ("y", -1)], "=", 0)
        r = h.solve("y", minimize=True)
        print(r)


    def skip_test3(self):
        """Attempt to adapt example from original code."""
        h = exact_simplex_helper.ExactSimplexHelper(["x", "y", "z"])
        h.add_constraint([("x", 15), ("y", 20), ("z", 25)], "<=", 1200)
        h.add_constraint([("x", 35), ("y", 60), ("z", 60)], "<=", 3000)
        h.add_constraint([("x", 20), ("y", 30), ("z", 25)], "<=", 1500)
        h.add_constraint([("y", 250)], ">=", 500)
        # FIXME the objective function isn't just one variable
        # c = list_to_fractions([300, 250, 450,  0,0,0,0])
        r = h.solve("y", minimize=False)
        print(r)

