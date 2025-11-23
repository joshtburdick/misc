# Test of flexible_lp solver.

import math
import unittest

import numpy as np

class Test(unittest.TestCase):

    def run1(self, m, n):
        """Runs a test case with m rows and n columns."""
        A = np.rand(   )



        self.assertLess( (x1 - x2).abs().mean() <= 1e-3 )

    def solve_scipy(self, ...):
        pass

    def solve_flexble_lp(self, ...):
        pass

    def test_lp_1(self):
        self.run1(3, 5)
        ...

if __name__ == '__main__':
    unittest.main()

