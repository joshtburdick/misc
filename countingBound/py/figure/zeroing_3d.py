#!/usr/bin/env python3
# Plots graphs, after zeroing out some hypercliques.
# FIXME
# - tweak axes to make this rectangular
# - make axes integers

import colorsys
import itertools
import math
import pdb

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import scipy.interpolate
from scipy.special import comb


# FIXME factor this out?
def scale_lightness(rgb, scale_l):
    """Scales the lightness of a color.

    From
    https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    """
    # convert RGB to HLS
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as RGB
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

class ZeroingPlot:
    """Plots effect of zeroing out one edge.

    """
    def __init__(self):
        # number of vertices
        self.n = 6
        # angle at which to place the 0'th vertex
        vertex_0_theta = 0
        all_cliques = list([frozenset(s) for s in itertools.combinations(range(self.n), 3)])
        # which edge will be zeroed out
        self.zeroed_edge = frozenset([0,1])
        # color cliques, depending on whether the edge hits it
        def color1(h):
            return colorsys.hsv_to_rgb(h, 0.5, 0.5)
        def color_clique(clique):
            zeroed_edge = self.zeroed_edge
            if zeroed_edge < clique:
                return color1(0)
            else:
                return color1(2/3)
        colors = {clique: color_clique(clique) for clique in all_cliques}
        # Set the vertices. Note that we angle these slightly, so that
        # they're "facing toward the origin". This is to emphasize
        # the importance of the "total number of cliques".
        theta = np.linspace(0, 2*np.pi, self.n, endpoint=False) + vertex_0_theta
        # useful constant
        c = np.sqrt(2) / 2
        self.vertex = np.stack([
            c * np.cos(theta), - c * np.cos(theta), np.sin(theta)])

    def get_location(self, cliques):
        """Gets the location of a clique."""
        # first, get the number "hit" by the zonked edge
        cliques_hit = [c for c in cliques if self.zeroed_edge in c]
        z = len(cliques_hit)
        n = len(cliques)
        return (n, z-n)

    def plot_num_functions(self):
        """Plots log_2(number of functions), as a surface."""
        log_num_functions = np.zeros([17, 5])
        for i in range(17):
            for j in range(5):
                # here, we're picking i cliques which aren't zonked,
                # and j cliques which are
                log_num_functions[i, j] = (math.log2(comb(16, i))
                    + math.log2(comb(4, j)))
        # plot this
        # ax = plt.axes(projection='3d')
        x = np.arange(5)
        y = np.arange(17)
        X, Y = np.meshgrid(x, y)
        Z = log_num_functions
        # pdb.set_trace()
        # fig = plt.figure()
        self.axs.plot_surface(X, Y, Z,
            rstride=1, cstride=1,
            cmap='binary', edgecolor='none', alpha=0.5)
        self.axs.set_xlabel('# cliques zonked')
        self.axs.set_ylabel('# cliques not zonked')
        self.axs.set_zlabel('lg(# functions)');

    def plot_clique_set(self, radius, center, cliques):
        """Plots a set of cliques.

        radius: the radius for the cliques
        z: the z-coordinate for the cliques
        cliques: the cliques, as a list of k-element tuples of ints
        Side effects: plots the cliques
        """
        # if there are no cliques, don't plot anything
        if not cliques:
            return
        # FIXME make this 3D
        for s in cliques:
            v = radius * self.vertex[:,list(s)] + np.array([center]).T
            plt.fill(v[0,:], v[1,:],
                edgecolor=self.colors[s],
                facecolor=scale_lightness(self.colors[s], 2),
                lw=3,
                alpha=self.alpha)

    def plot_clique_sets(self):
        pass

    def plot_it(self):
        """Plots the rectangle containing all the functions."""
        self.fig = plt.figure(figsize=(9, 5))
        self.axs = self.fig.add_subplot(111, projection='3d')
        # axs = plt.axes()
        # plot (log_2 of) the number of functions
        self.plot_num_functions()

        # plot some random sets of cliques
        cliques1 = frozenset([frozenset([0,1,3]), frozenset([2,3,4])])
        # self.plot_clique_set(0.5, 0, cliques1)

        # for practice: draw a triangle
        # vertices = np.array([[0,0,0], [1,0,0], [0,1,0]])
        vertices = [list(zip([0,0,0],[4,0,0],[0,4,0]))]
        poly = Poly3DCollection(vertices, alpha=0.8, color='green')
        self.axs.add_collection3d(poly)

        plt.savefig('zeroing_3d.png')  # , bbox_inches='tight')


if __name__ == '__main__':
    z = ZeroingPlot()
    z.plot_it()

