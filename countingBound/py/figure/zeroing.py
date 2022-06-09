#!/usr/bin/env python3
# Plots graphs, after zeroing out some hypercliques.

import colorsys
import itertools
import pdb

import matplotlib
import matplotlib.pyplot as plt
import numpy as np



def scale_lightness(rgb, scale_l):
    """Scales the lightness of a color.

    From
    https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    """
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

class CliqueFigure:
    """Plots cliques in a figure.

    This assumes that a given hyperedge is the same color
    in each plot.
    """

    def __init__(self, n, colors, vertex_0_theta):
        """
        n: number of vertices in each graph
        colors: hash from sets (defined as sorted lists of k numbers)
            to color of each set
        vertex_0_theta: angle at which to place vertex 0
        """
        self.n = n
        self.colors = colors
        self.alpha = 0.5
        # set the vertices
        theta = np.linspace(0, 2*np.pi, n, endpoint=False) + vertex_0_theta
        self.vertex = np.stack([np.cos(theta), np.sin(theta)])

    def plot_sets(self, radius, center, sets):
        """Plots some sets.

        radius: the radius for the sets
        center: the center, as a 2-element tuple, or np.array of shape (2,)
        sets: the sets, as a list of k-element tuples of ints
        Side effects: plots the sets
        """
        for s in sets:
            v = radius * self.vertex[:,s] + np.array([center]).T
            plt.fill(v[0,:], v[1,:],
                edgecolor=self.colors[s],
                facecolor=scale_lightness(self.colors[s], 2),
                lw=3,
                alpha=self.alpha)





def plot_Z_relation():
    """Plots the 'zeroing-one-edge' relation."""
    plt.figure(figsize=(6,7))
    plt.xlim(-1,7)
    plt.ylim(-1, 5)
    def color1(h):
        return colorsys.hsv_to_rgb(h, 0.5, 0.5)
    colors = {
        (0,1,2):color1(0/4),
        (0,1,3):color1(1/4),
        (0,2,3):color1(2/4),
        (1,2,3):color1(3/4)}
    cf = CliqueFigure(4, colors, 0)
    # lay out coordinates for each set; this will be keyed by set,
    # and its value will be coordinates
    set_location = {} 
    edges = list(itertools.combinations(range(4), 3))
    for i in range(5):
        print(i)
        subsets = tuple(itertools.combinations(edges, i))
        for j in range(len(subsets)):
            print(((i,j), subsets[j]))
            set_location[subsets[j]] = (i,j)
    # plot effects of zeroing out an edge



    # cf.plot_sets(0.4, np.array([0,0.1]), [(0,1,2), (0,1,3)])
    plt.savefig('Z.png')

def plot_zeroing():
    """Plots effect of zeroing out one edge."""
    pass









if __name__ == '__main__':
    plot_Z_relation()

