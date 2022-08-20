#!/usr/bin/env python3
# Plots graphs, after zeroing out some hypercliques.
# FIXME
# - show which edge was zeroed!
# - tweak alpha for cliques which are hit by an edge?

import colorsys
import itertools
import pdb
import random

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

def scale_lightness(rgb, scale_l):
    """Scales the lightness of a color.

    From
    https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    """
    # convert RGB to HLS
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as RGB
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
        FIXME
        - include radius here?
        - add method to show which edge was zeroed out?
        """
        self.n = n
        self.colors = colors
        self.alpha = 0.4
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
            v = radius * self.vertex[:,list(s)] + np.array([center]).T
            plt.fill(v[0,:], v[1,:],
                edgecolor=self.colors[s],
                facecolor=scale_lightness(self.colors[s], 2),
                lw=1,
                alpha=self.alpha,
                joinstyle='round')

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
        np.mean(endpoints, axis=0) + (0, 0.4),
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
        # note that < is comparing frozensets
        return frozenset([c for c in cliques if not frozenset(edge) < c])
    # note that this removes self-loops
    return list(set([cliques_remaining(e) for e in edges]) - set([cliques]))

def plot_Z_relation():
    """Plots the 'zeroing-one-edge' relation."""
    plt.figure(figsize=(6,7))
    plt.axis('off')
    plt.xlim(-3, 3.5)
    plt.ylim(0, 4.5)
    def color1(h):
        return colorsys.hsv_to_rgb(h, 0.5, 0.5)
    colors = {
        frozenset((0,1,2)):color1(0/4),
        frozenset((0,1,3)):color1(1/4),
        frozenset((0,2,3)):color1(2/4),
        frozenset((1,2,3)):color1(3/4)}
    cf = CliqueFigure(4, colors, 0)
    # lay out coordinates for each set; this will be keyed by set,
    # and its value will be coordinates
    set_location = {} 
    all_cliques = list([frozenset(s) for s in itertools.combinations(range(4), 3)])
    for j in range(0, 5):
        print(j)
        # ??? should this be a set rather than a tuple?
        subsets = tuple(itertools.combinations(all_cliques, j))
        for i in range(len(subsets)):
            print(((i,j), subsets[i]))
            # this is mostly centered, but also slightly tilted
            set_location[frozenset(subsets[i])] = (i - len(subsets)/2 + j/3, j)
    # plot effects of zeroing out an edge
    for (cliques, location) in set_location.items():
        cliques_below = zero_out_edges(cliques)
        for c in cliques_below:
            location_1 = set_location[c]
            p = curved_line(np.array([location, location_1]))
            plt.plot(p[:,0], p[:,1], c='black', alpha=0.3)

    # plot the sets
    for (cliques, location) in set_location.items():
        # pdb.set_trace()
        print("cliques =")
        print(cliques)
        print("location = " + str(location))
        cf.plot_cliques(0.25, location, cliques)
    # cf.plot_sets(0.4, np.array([0,0.1]), [(0,1,2), (0,1,3)])
    plt.savefig('Z.pdf', bbox_inches='tight')

def plot_zeroing_rectangle():
    """Plots rectangle showing effect of zeroing out one edge."""
    n = 6
    # which edge will be zeroed out
    zeroed_edge = frozenset([0,1])
    plt.figure(figsize=(8,3))
    # plt.axis('off')
    margin = 0.5
    plt.xlim(-margin, 16+margin)
    plt.ylim(-margin, 4+margin)
    plt.xlabel('# cliques not zonked')
    plt.ylabel('# cliques zonked')
    # list of all cliques
    all_cliques = [frozenset(s) for s in itertools.combinations(range(n), 3)]
    # which cliques are hit, and which are missed
    hit_cliques = [c for c in all_cliques if zeroed_edge < c]
    missed_cliques = [c for c in all_cliques if not zeroed_edge < c]

    # color scheme: red if the clique is hit, otherwise blue
    def color1(h):
        return colorsys.hsv_to_rgb(h, 0.5, 0.5)
    colors = {c: colorsys.hsv_to_rgb(0, 1, 0.5) if zeroed_edge < c
            else colorsys.hsv_to_rgb(2/3, 0.5, 0.5)
        for c in all_cliques}
    cf = CliqueFigure(6, colors, 0)

    # pick a random set of cliques "representative" of each point in the rectangle
    for x in range(17):
        for y in range(5):
            cliques = random.sample(hit_cliques, y) + random.sample(missed_cliques, x)
            cf.plot_cliques(0.3, (x,y), cliques)
            # also, draw a line

    # add lines showing "levels" of set size
    # XXX this could be done with clipping, but this hopefully will also work,
    # and is fairly easy to think about
    for x in range(1, 17):
        for y in range(0, 4):
            plt.plot([x, x-1], [y, y+1], c='black', lw=0.5, alpha=0.5)

    plt.savefig('zeroingRectangle.pdf', bbox_inches='tight')

if __name__ == '__main__':
    plot_Z_relation()
    plot_zeroing_rectangle()

