#!/usr/bin/env python3
# Plots graphs, after zeroing out some hypercliques.

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
        center: the center, as a 2-element numpy array
        sets: the sets, as a list of k-element tuples of ints
        Side effects: plots the sets
        """
        for s in sets:
            v = radius * self.vertex[:,s] + np.array([center]).T
            plt.fill(v[0,:], v[1,:], facecolor=self.colors[s], alpha=self.alpha)






def plot_it():
    plt.figure()
    colors = {(0,1,2):'red', (0,1,3):'orange', (0,2,3):'green', (1,2,3):'blue'}
    cf = CliqueFigure(4, colors, 0)
    cf.plot_sets(1, np.array([0,0.1]), [(0,1,2), (0,2,3), (1,2,3)])
    plt.savefig('test0.png')

if __name__ == '__main__':
    plot_it()

