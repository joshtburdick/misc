#!/usr/bin/env python3
# Plots graphs, after zeroing out some hypercliques.

import colorsys
import itertools
import math
import pdb

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
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
        # ??? This has some hard-coded constants. But the number of
        # vertices may not change much.
        self.n = 6
        # angle at which to place the 0'th vertex
        vertex_0_theta = 0
        # radius of the cliques
        self.radius = 0.6
        self.alpha = 0.5
        self.all_cliques = list([frozenset(s) for s in itertools.combinations(range(self.n), 3)])
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
        self.colors = {clique: color_clique(clique)
            for clique in self.all_cliques}
        # Set the vertices. Note that we angle these slightly, so that
        # they're "facing toward the origin". This is to emphasize
        # the importance of the "total number of cliques".
        theta = np.linspace(0, 2*np.pi, self.n, endpoint=False) + vertex_0_theta
        # useful constant
        c = np.sqrt(2) / 2
        self.vertex = np.stack([
            c * np.cos(theta), - c * np.cos(theta), np.sin(theta)], axis=1)

    def get_location(self, cliques):
        """Gets the location of a clique."""
        # first, get the number "hit" by the zonked edge
        cliques_hit = [c for c in cliques if self.zeroed_edge < c]
        z = len(cliques_hit)
        n = len(cliques)
        return (z, n-z)

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
            cmap='binary', edgecolor='none', alpha=0.3)
        self.axs.set_xlabel('# cliques zonked')
        self.axs.set_ylabel('# cliques not zonked')
        self.axs.set_zlabel('lg(# functions)');

    def plot_clique_set(self, center, cliques):
        """Plots a set of cliques.

        center: the center of the cliques
        cliques: the cliques, as a list of k-element tuples of ints
        Side effects: plots the cliques
        """
        # if there are no cliques, don't plot anything
        if not cliques:
            return
        # loop through the cliques
        print(center)
        for c in cliques:
            i = [x for x in c]
            v = self.radius * self.vertex[i,:] + center
            vertices1 = [[tuple(r) for r in v]]
            poly = Poly3DCollection(vertices1, alpha=self.alpha,
                color=self.colors[c])
            self.axs.add_collection3d(poly)
            # plt.fill(v[0,:], v[1,:], v[2,:],
            #     edgecolor=self.colors[c],
            #     facecolor=scale_lightness(self.colors[c], 2),
            #    lw=3,
            #    alpha=self.alpha)

    def plot_clique_sets(self, clique_sets):
        stacking_height = 2
        # this is the height of the next clique to be plotted at
        # particular coordinates (so far)
        z = np.zeros([5, 17])
        # loop through the sets of cliques
        for s in clique_sets:
            (x, y) = self.get_location(s)
            # put this clique above any other cliques in this column
            center = (x, y, z[x,y])
            z[x,y] += stacking_height
            self.plot_clique_set(center, s)

        # lastly, draw a line connecting all the sets of cliques
        # in this "stack"?
        # FIXME

    def plot_it(self):
        """Plots the rectangle containing all the functions."""
        # set up axes
        self.fig = plt.figure(figsize=(9, 5))
        self.axs = self.fig.add_subplot(111, projection='3d')
        self.axs.set_xlim(4,0)
        self.axs.xaxis.set_major_locator(MaxNLocator(integer=True))
        # make x- and y-axis scales be square, and flatten it somewhat
        self.axs.set_box_aspect((2,4,3))
        # axs = plt.axes()
        # plot (log_2 of) the number of functions
        self.plot_num_functions()

        # plot some random sets of cliques
        cliques1 = [frozenset([frozenset([0,1,3]), frozenset([2,3,4])])]
        cliques1.append(self.all_cliques)
        for i in range(5):
            cliques1.append(frozenset(np.random.choice(self.all_cliques, 1)))
            cliques1.append(frozenset(np.random.choice(self.all_cliques, 3)))
            cliques1.append(frozenset(np.random.choice(self.all_cliques, 5)))
            cliques1.append(frozenset(np.random.choice(self.all_cliques, 10)))
            cliques1.append(frozenset(np.random.choice(self.all_cliques, 15)))

        self.plot_clique_sets(cliques1)

        # for practice: draw a triangle
        # vertices0 = np.array([[0,0,0], [4,0,0], [0,4,0]])
        # vertices = [list(zip([0,0,0],[4,0,0],[0,4,0]))]
        # pdb.set_trace()
        # poly = Poly3DCollection(vertices, alpha=0.5, color='green')
        # self.axs.add_collection3d(poly)

        plt.savefig('zeroing_3d.png')  # , bbox_inches='tight')


if __name__ == '__main__':
    z = ZeroingPlot()
    z.plot_it()

