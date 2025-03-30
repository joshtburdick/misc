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

# plt.rcParams["text.usetex"] = True

class BouncePlot:

    def __init__(self, n, k):
        self.n = n
        self.k = k
        # number of possible cliques
        self.N = scipy.special.comb(n, k)
        # number of cliques possibly "hit" by zeroing out an edge
        self.h = scipy.special.comb(n-2, k-2)
        self.rng = np.random.default_rng(12345)

    def sample_step(self, U_t0, num_samples=1000):
        """Samples one 'bounce'.

        U_t0: starting point (x, y) = (num. cliques, num. gates)
        Returns: a dict with keys "U_t0", "V_t1", and "U_t1".
            U_t2 and V_t2 will be numpy arrays with columns "x", "y_lo", and "y_hi".
            ... also, "V_t1_bounds" and "U_t2_bounds", which are endpoints of
                polygons indicating the bounds.
        """
        U_t0 = np.array(U_t0)
        # bounce "down"
        num_hit = self.rng.hypergeometric(U_t0[0], self.N - U_t0[0], self.h, num_samples)
        V_t1_x = U_t0[0] - num_hit
        V_t1_y_lo = U_t0[1] - num_hit
        V_t1_y_hi = np.full(num_samples, [U_t0[1]])

        max_h = min(U_t0[0], self.h)
        # note that we duplicate the last point, to streamline drawing connecting lines
        # ??? maybe not needed
        V_t1_bounds = np.stack([U_t0, U_t0 - (max_h, max_h), U_t0 - (max_h, 0)])

        # bounce "up"
        num_added = self.rng.binomial(self.h, 0.5, num_samples)
        U_t2_x = V_t1_x + num_added
        U_t2_y_lo = V_t1_y_lo
        U_t2_y_hi = V_t1_y_hi + num_added
        U_t2_max_x = min(self.N, U_t0[0]+self.h)
        # amount which this could go "to the right" of the starting point
        # (which is bounded by the left and right edges of the plot)
        dx = min(min(self.h, self.N - U_t0[0]), U_t0[0])
        # XXX this is finicky, and possibly incorrect
        U_t2_bounds = np.stack([
            V_t1_bounds[1,:],
            V_t1_bounds[2,:],
            V_t1_bounds[2,:] + (self.h, self.h),
            V_t1_bounds[2,:] + (self.h + dx, self.h),
            V_t1_bounds[1,:] + (self.h + dx, dx),
            V_t1_bounds[1,:] + (self.h, 0)])

        return {
            "U_t0": U_t0,
            "V_t1": np.stack([V_t1_x, V_t1_y_lo, V_t1_y_hi], axis=1),
            "U_t2": np.stack([U_t2_x, U_t2_y_lo, U_t2_y_hi], axis=1),
            "V_t1_bounds": V_t1_bounds,
            "U_t2_bounds": U_t2_bounds
        }

    def plot_lines(self, axs, p, hue, alpha=0.4):
        """Plot lines.

        FIXME
        - show lines between sampling from a given starting point?

        axs: Axes object to draw on
        p: endpoints, as a numpy array with three columns
        hue: the hue
        """
        rgb = matplotlib.colors.hsv_to_rgb(np.array([hue, 1, 0.75]))

        # reshape points
        p1 = [ [(p[i,0], p[i,1]), (p[i,0], p[i,2])]
            for i in range(p.shape[0]) ]
        # pdb.set_trace()
 
        lines = matplotlib.collections.LineCollection(
             p1, colors=np.concat([rgb, np.array([alpha])]))
        axs.add_collection(lines)


    def plot_bounce(self, axs, x_proportion, num_samples=100):
        """Plots starting with some number of cliques."""
        # XXX connect lines showing trajectories of individual samples?
        x = int(self.N * x_proportion)
        s = self.sample_step((x, 0), num_samples=num_samples)
        # pdb.set_trace()
        axs.scatter(x, 0, color="black", s=8, alpha=1)
        alpha = 0.02
        self.plot_lines(axs, s["V_t1"], 0, alpha=alpha)
        self.plot_lines(axs, s["U_t2"], 2/3, alpha=alpha)

        # plot possible regions
        axs.fill(
            s["V_t1_bounds"][:,0],
            s["V_t1_bounds"][:,1],
            facecolor="#ff000020", edgecolor="#ff0000ff", linewidth=0.025)
        axs.fill(
            s["U_t2_bounds"][:,0],
            s["U_t2_bounds"][:,1],
            facecolor="#0000ff20", edgecolor="#0000ffff", linewidth=0.025)
        axs.set_xlabel("Number of cliques")
        axs.set_ylabel("Relative number of gates")

# FIXME force this to have 1:1 aspect ratio?
plt.figure(figsize=(8,2.5))
# bp = BouncePlot(12, 5)
bp = BouncePlot(14, 4)

for p in [0, 0.25, 0.5, 0.75, 1]:
    bp.plot_bounce(plt.gca(), p)

plt.margins(0.02)
plt.tight_layout()
plt.savefig("../../bound2/walking_bounds_0.pdf")

