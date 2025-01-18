#!/usr/bin/env python3
# Plots bounds implied by random walk.


import colorsys
import itertools
import pdb
import random

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import scipy.interpolate
import scipy.special
import scipy.stats

class BouncePlot:

    def __init__(self, n, k):
        self.n = n
        self.k = k
        # number of possible cliques
        self.N = scipy.special.comb(n, k)
        # number of cliques "hit" by zeroing out an edge
        self.h = scipy.special.comb(n-2, k-2)
        self.rng = np.random.default_rng(12345)

    def sample_step(self, U_t0, num_samples=1000):
        """Samples one 'bounce'.

        U_t0: starting point (x, y) = (num. cliques, num. gates)
        Returns: a dict with keys "U_t0", "V_t1", and "U_t1".
            U_t2 and V_t2 will be numpy arrays with columns "x", "y_lo", and "y_hi".
        """
        # bounce "down"
        num_hit = self.rng.hypergeometric(U_t0[0], self.N - U_t0[0], self.h, num_samples)
        V_t1_x = U_t0[0] - num_hit
        V_t1_y_lo = U_t0[1] - num_hit
        V_t1_y_hi = np.full(num_samples, [U_t0[1]])
        
        # bounce "up"
        num_added = self.rng.binomial(self.h, 0.5, num_samples)
        U_t2_x = V_t1_x + num_added
        U_t2_y_lo = V_t1_y_lo + num_added
        U_t2_y_hi = V_t1_y_hi + num_added

        return {
            "U_t0": U_t0,
            "V_t1": np.stack([V_t1_x, V_t1_y_lo, V_t1_y_hi]),
            "U_t2": np.stack([U_t2_x, U_t2_y_lo, U_t2_y_hi]),
        }

bp = BouncePlot(10, 3)

z = bp.sample_step((30, 0), num_samples=3)

pdb.set_trace()

