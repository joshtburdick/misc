#!/usr/bin/env python3
# Plots graphs, after zeroing out some hypercliques.

import colorsys
import itertools
import pdb

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate


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

    def plot_cliques(self, radius, center, cliques):
        """Plots some cliques.

        radius: the radius for the cliques
        center: the center, as a 2-element tuple, or np.array of shape (2,)
        cliques: the cliques, as a list of k-element tuples of ints
        Side effects: plots the cliques
        """
        # if there are no cliques, don't plot anything
        if not cliques:
            return
        for s in cliques:
            v = radius * self.vertex[:,s] + np.array([center]).T
            plt.fill(v[0,:], v[1,:],
                edgecolor=self.colors[s],
                facecolor=scale_lightness(self.colors[s], 2),
                lw=3,
                alpha=self.alpha)

def interpolate(nodes, num_points=101):
    """Does spline interpretation.

    Essentially taken from
    https://stackoverflow.com/questions/29837854/matplotlib-draw-spline-from-multiple-points
    nodes: points to interpolate between (as a numpy array,
        with one row per point)
    num_points: number of interpolated points to include
    Returns: a numpy array with that many rows
    """
    x = nodes[:,0]
    y = nodes[:,1]
    tck,u = scipy.interpolate.splprep([x,y], k=2)
    xnew, ynew = scipy.interpolate.splev(np.linspace(0, 1,num_points), tck,der = 0)
    return np.stack([xnew, ynew], axis=1)

def curved_line(endpoints):
    """Draws a curved line between the endpoints."""
    p = np.array([
        endpoints[0,:],
        # in the middle is the midpoint, except a bit higher up
        np.mean(endpoints, axis=0) + (0, 0.8),
        endpoints[1,:] ])
    return interpolate(p)

def zero_out_edges(cliques):
    """
    Finds the cliques remaining, after zeroing out an edge.

    """
    # get the vertices which are relevant
    vertices = set(list(itertools.chain(*cliques)))
    # the edges, which we could zero out
    edges = itertools.combinations(vertices, 2)
    def cliques_remaining(edge):
        return tuple([c for c in cliques if not set(edge) < set(c)])
    return list(tuple([cliques_remaining(e) for e in edges]))

def plot_Z_relation():
    """Plots the 'zeroing-one-edge' relation."""
    plt.figure(figsize=(6,7))
    plt.xlim(-3.5, 3.5)
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
    all_edges = list(itertools.combinations(range(4), 3))
    for j in range(0, 5):
        print(j)
        # ??? should this be a set rather than a tuple?
        subsets = tuple(itertools.combinations(all_edges, j))
        for i in range(len(subsets)):
            print(((i,j), subsets[i]))
            set_location[subsets[i]] = (i + 1/2 - len(subsets)/2, j)
    # plot effects of zeroing out an edge
    for (cliques, location) in set_location.items():
        cliques_below = zero_out_edges(cliques)
        for c in cliques_below:
            location_1 = set_location[c]
            p = curved_line(np.array([location, location_1]))
            plt.plot(p[:,0], p[:,1], c='black', alpha=0.3)

 
    # plot the sets
    for (cliques, location) in set_location.items():
        cf.plot_cliques(0.25, location, cliques)
    # cf.plot_sets(0.4, np.array([0,0.1]), [(0,1,2), (0,1,3)])
    plt.savefig('Z.png')

def plot_zeroing():
    """Plots effect of zeroing out one edge."""
    pass








if __name__ == '__main__':
    plot_Z_relation()

