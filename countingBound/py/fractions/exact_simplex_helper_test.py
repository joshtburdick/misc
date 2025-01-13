# Some small tests.
# Those starting with _test currently don't work (and so are skipped).

import unittest

import fractions

import exact_simplex_helper

class TestExactSimplexHelper(unittest.TestCase):


    def _test1(self):
        """Example from explication."""
        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2"])
        h.add_constraint([("x1", 1), ("x2", 2)], "<=", 4)
        h.add_constraint([("x1", 1), ("x2", -1)], "<=", 1)
        r = h.solve_1([("x1", 3), ("x2", 2)])
        self.assertEqual(r["x1"], 2)
        self.assertEqual(r["x2"], 1)
        print(r)


    def test2(self):
        """??? not sure why this one is failing? """
        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2", "y"])
        h.add_constraint([("x1", 1), ("x2", 5)], "<=", 10)
        h.add_constraint([("x1", 2), ("x2", 1)], "<=", 10)
        h.add_constraint([("x1", 1), ("x2", 1), ("y", -1)], "=", 0)
        r = h.solve("y", minimize=False)
        print(r)


    def _test3(self):
        """Attempt to adapt example from original code."""
        h = exact_simplex_helper.ExactSimplexHelper(["x", "y", "z"])
        h.add_constraint([("x", 15), ("y", 20), ("z", 25)], "<=", 1200)
        h.add_constraint([("x", 35), ("y", 60), ("z", 60)], "<=", 3000)
        h.add_constraint([("x", 20), ("y", 30), ("z", 25)], "<=", 1500)
        h.add_constraint([("y", 250)], ">=", 500)
        r = h.solve_1([("x", 300), ("y", 250), ("z", 450)])
        print(r)

    def _test_from_wiki_1_tweaked(self):
        """Modified example from Wikipedia.

        I had been thinking that the initial solution had to have
        all variables >= 0. This doesn't seem to always be the case.
        """
        h = exact_simplex_helper.ExactSimplexHelper(["x", "y", "z"])
        h.add_constraint([("x", 3), ("y", 2), ("z", 1)], "<=", 10)
        h.add_constraint([("x", 2), ("y", 5), ("z", 3)], "<=", 15)
        h.add_constraint([("x", 1)], ">=", 1)
        h.add_constraint([("y", 1)], ">=", 1)
        h.add_constraint([("z", 1)], ">=", 1)
        # we're minimizing this, so we negate these coefficients
        r = h.solve_1([("x", 2), ("y", 3), ("z", 4)])
        print(r)


    def test_crossed_lines(self):
        """Test based on multiple lines.

        This prints x = y = 0 as a solution, which is incorrect.
        But this may reproduece the previously-seen problem...
        """
        print("-------- start crossed lines test")
        h = exact_simplex_helper.ExactSimplexHelper(["x", "y"])
        h.add_constraint([("x", 4), ("y", 1)], ">=", 4)
        h.add_constraint([("x", fractions.Fraction(3,2)), ("y", 2)], ">=", 3)
        h.add_constraint([("x", fractions.Fraction(2,3)), ("y", 3)], ">=", 2)
        h.add_constraint([("x", fractions.Fraction(1,4)), ("y", 4)], ">=", 1)
        r = h.solve_1([("x", -1), ("y", -1)])
        print(r)
        print("FIXME check if solution is correct")
        print("-------- end crossed lines test")


    def test_from_wiki_2(self):
        """Example from Wikipedia which possibly needs two phases."""
        print("in test_from_wiki_2()...")
        h = exact_simplex_helper.ExactSimplexHelper(["x", "y", "z"])
        h.verbosity = 3
        h.add_constraint([("x", 3), ("y", 2), ("z", 1)], "=", 10)
        h.add_constraint([("x", 2), ("y", 5), ("z", 3)], "=", 15)
        # we're minimizing this, so we negate these coefficients
        r = h.solve_1([("x", 2), ("y", 3), ("z", 4)])
        self.assertEqual(r["x"], fractions.Fraction(15,7))
        self.assertEqual(r["y"], 0)
        self.assertEqual(r["z"], fractions.Fraction(25,7))


    def test_two_phase(self):
        """
        Two-phase example from
        https://sites.math.washington.edu/~burke/crs/407/notes/section3-18.pdf
        """
        h = exact_simplex_helper.ExactSimplexHelper(["x1", "x2", "x3"])
        h.add_constraint([("x1", 2), ("x2", -1), ("x3", 2)], "<=", 4)
        h.add_constraint([("x1", 2), ("x2", -3), ("x3", 1)], "<=", -5)
        h.add_constraint([("x1", -1), ("x2", 1), ("x3", -2)], "<=", -1)
        # here, we're maximizing
        r = h.solve_1([("x1", 1), ("x2", -1), ("x3", 1)])
        print(r)
        self.assertEqual(r["x1"], 0)
        self.assertEqual(r["x2"], fractions.Fraction(14,5))
        self.assertEqual(r["x3"], fractions.Fraction(17,5))


    def test_silly_1(self):
        """A really silly test."""
        print("in test_silly_1:")
        h = exact_simplex_helper.ExactSimplexHelper(["x"])
        h.add_constraint([("x", 1)], ">=", 1)
        # we're minimizing this, so we negate these coefficients
        r = h.solve_1([("x", -1)])
        print(r)
        # what the result should be
        self.assertEqual(r["x"], 1)


    def test_silly_2(self):
        """Another really silly test."""
        print("in test_silly_2:")
        h = exact_simplex_helper.ExactSimplexHelper(["x"])
        h.add_constraint([("x", 1)], "=", 1)
        # we're minimizing this, so we negate these coefficients
        r = h.solve_1([("x", -1)])
        print(r)
        # what the result should be
        self.assertEqual(r["x"], 1)


