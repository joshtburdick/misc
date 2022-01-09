#!/usr/bin/env python3
# Estimates rank of finding various numbers of functions.
# These try to emphasize clarity over efficiency.

import matplotlib
import matplotlib.pyplot as plt
import pdb

import numpy as np
import scipy.special
from scipy.special import binom, comb
import scipy.stats
from scipy.stats import hypergeom

def rankBoundZeroedEdges(n, k):
    """Estimate of rank for finding some number of cliques.

    This version "zonks" edges.
    n: size of the input graph
    k: size of the cliques to find
    Returns: a vector r of length N+1, where N = choose(n, k),
        such that r_i is an estimate of the rank (in the list of all
        circuits) of the smallest circuit finding i cliques.
    """
    # number of cliques
    numCliques = comb(n, k, exact=True)
    # the bound (initially all 0)
    r = np.zeros([numCliques+1])
    # probability of a clique being "missed" when you feed in a 0
    # to a given edge
    p = 1 - (comb(k,2) / binom(n,2))
    for i in range(1, numCliques+1):
        # This "weight" is the distribution of smaller sets of cliques
        # which the current level definitely requires more gates than.
        # Note that the 0 can "miss" all of the cliques; in this case,
        # our lower bound is 0 (which is what r is initialized to).
        # (As i grows, this case will be increasingly rare.)
        w = scipy.stats.binom.pmf(range(i), i, p)
        # The bound, at the i'th level, is the expected size of the
        # circuits which are smaller, plus half the number of new
        # circuits we're adding at this level.
        r[i] = sum(w * r[0:i]) + comb(numCliques, i) / 2
    return r

def rankBoundZeroedVertices(n, k):
    """Estimate of rank for finding some number of cliques.

    This version feeds in zeros to all edges incident to a vertex.
    n: size of the input graph
    k: size of the cliques to find
    Returns: a dict r, with keys (numVertices, numCliques),
        giving a lower bound on the expected rank of the smallest circuit
        finding numCliques cliques in a graph with <= numVertices vertices.
    """
    # number of cliques
    # numCliques = scipy.special.comb(n, k, exact=True)
    # the bound: initially there are only two functions
    # FIXME make this a dict of numpy vectors, or a ragged array?
    r = {(k,0):0, (k,1):1}
    # loop through the number of vertices
    for v in range(k+1, n+1):
        # loop through the number of cliques with that many vertices
        maxCliques = comb(v, k, exact=True)
        for j in range(maxCliques+1):
            # the number of cliques which include an arbitrary vertex
            numCliquesIncludingVertex = comb(v-1, k-1, exact=True)
            # the number of cliques potentially zonked by
            # zeroing out edges connected to a vertex
            maxCliquesZonked = min(j, numCliquesIncludingVertex)
            # The number of cliques zonked. Note that we require at least
            # one clique to be zonked; this may not be necessary.
            z = range(1, maxCliquesZonked+1)
            # the probability of some number of these being zonked
            w = hypergeom(
                    # number of possible cliques
                    maxCliques,
                    # number of those present
                    j,
                    # number of possible cliques we're sampling
                    numCliquesIncludingVertex
                    ).pmf(z)
            # the "expected value" bound for this many cliques
            # (defined on this many vertices)
            r[(v,j)] = (
                    # the expected rank of functions which are definitely smaller
                    # (a weighted sum of the number of cliques "left over" after
                    # zonking cliques from one vertex) ...
                    sum( w * np.array([r[(v,j-z1)] for z1 in z]) )
                    # ... plus half the number of functions (with this many
                    # cliques, defined on this many vertices) being added
                    # ??? should w also be weighting this ? I think that it
                    # doesn't matter, by linearity of expectation.
                    + comb(maxCliques, j) / 2)
    return r

def rankBoundZeroingVertexEdges(maxNumVertices, k):
    """Bounds function rank, zeroing a vertex' edges, one at a time.

    This version adds vertices, one at a time. For each, it adds
    edges one at a time, and tracks how many cliques would be lost
    by zeroing out that edge.
    XXX this function is getting looooong

    maxNumVertices: size of the input graph
    k: size of the cliques to find
    Returns: a dict r, keyed by (v, e), where:
        v: the number of vertices "completely" added so far
            (which ranges up to maxNumVertices)
        e: the number of edges added to the "new" vertex
            (this can be 0)
        The value r[(v,e)] is a numpy array such that r[(v,e)][i] is a lower
        bound on E[rank] of all functions with exactly i cliques.
    """
    # the bound: initially there are only two functions, defined on
    # k vertices (either there's a k-clique, or not)
    # ??? rename this?
    # ??? should the rank in the bound be 0-based? does it matter?
    bound = { (k,0): np.array([0, 1]) }
    # loop through the number of vertices, including the new vertex v
    for v in range(k+1, maxNumVertices+1):
        # the bound, without the new vertex (initially with no edges added)
        bound1 = bound[(v-1, 0)]
        # Loop through the number of edges connected to v.
        # Here, j is the number of edges after adding a new edge, e.
        # (Note that when we add a new vertex, we initially need
        # at least k-1 edges, in order for the new vertex to possibly
        # be part of a k-clique).
        for numNewEdges in range(k-1, v-1):
            # number of cliques so far (formed from vertices other
            # than v, plus edges connected from v so far)
            numCurrentCliques = bound1.shape[0] - 1
            # The number of cliques created by adding an edge e;
            # all of these would be "zonked" if e were set to 0.
            # e's ends determine two vertices of these cliques, but
            # the remaining edges could be any subset of the remaining edges.
            numNewCliques = comb(numNewEdges-1, k-2, exact=True)
            # the total number of cliques which are possible
            # (with the new edge added)
            maxCliques = numCurrentCliques + numNewCliques
            # the bound on rank for every possible number of cliques
            # (with the new edge, more cliques are possible)
            bound2 = np.zeros(maxCliques+1)
            # numCliques is the total number of cliques in the graph,
            # with the edge added (whether they include it, or not)
            for numCliques in range(0, maxCliques+1):
                # the number of cliques potentially zonked by
                # zeroing out edges connected to a vertex
                maxCliquesZonked = min(numCliques, numNewCliques)
                # The number of cliques zonked. Note that we don't
                # require any cliques to be zonked (this is basically
                # including the functions from the previous step).
                z = np.arange(maxCliquesZonked+1)
                # the probability of some number of these being zonked
                w = hypergeom(
                    # number of possible cliques
                    maxCliques,
                    # number of those present
                    numCliques,
                    # number of cliques which intersect the new edge
                    # (and so could be zonked)
                    numNewCliques
                    ).pmf(z)
                # The number of _new_ functions with this many cliques.
                # This is the total with the new edge added, minus
                # the total without the new edge (this latter number
                # can be zero).
                numNewFunctions = (comb(maxCliques, numCliques, exact=True)
                    - comb(numCurrentCliques, numCliques, exact=True))
                # the bound on E[rank] for this many cliques
                bound2[numCliques] = (
                        # first, the bound from the "previous" cliques,
                        # with the edge zonked
                        w.dot(bound1[numCurrentCliques-z])
                        # plus, half the new functions 
                        + numNewFunctions / 2)
            # update the bound
            bound1 = bound2
            # FIXME also save this when the new vertex is only
            # partially added (mostly for debugging) ?
        bound[(v,0)] = bound2
    return bound

def plotBoundAtLevels(n, k):
    """Plots the bound at various levels.

    n: number of vertices
    k: number of vertices in cliques
    Side effects: plots different bounds on the rank of the function
        finding some number of cliques.
    """
    print('n = ' + str(n) + ', k = ' + str(k))
    maxCliques = comb(n, k, exact=True)
    # compute the bounds
    bound1 = rankBoundZeroedEdges(n, k)
    bound2All = rankBoundZeroedVertices(n, k)
    bound2 = [bound2All[n,k] for k in range(maxCliques+1)]
    # plot
    plt.figure()
    plt.plot(range(maxCliques+1), bound1, label='Zeroing out edges')
    # ??? should this be included?
    plt.plot(range(maxCliques+1), bound2, label='Zonking vertices')
    plt.plot(range(maxCliques+1),
            [comb(maxCliques, c)/2 for c in range(maxCliques+1)],
            label='Naive counting bound')
    plt.title('n = ' + str(n) + ', k = ' + str(k))
    plt.xlabel('Number of cliques')
    plt.ylabel('E[rank of functions]')
    plt.legend()
    # force x-axis to be plotted as integers
    plt.gca().xaxis.get_major_locator().set_params(integer=True)
    # plt.savefig('rank_n=' + str(n) + '_k=' + str(k) + '.pdf')
    plt.savefig('rank_n=' + str(n) + '_k=' + str(k) + '.png')

if False:
    for n in range(6, 12):
        plotBoundAtLevels(n, 3)
        plotBoundAtLevels(n, 4)
        plotBoundAtLevels(n, 5)

for k in range(3, 7):
    plotBoundAtLevels(2*k, k)

# this runs, but gives a not-great bound
# z = rankBoundZeroingVertexEdges(6,3)
# pdb.set_trace()

if False:
    plotBoundAtLevels(6, 3)
    plotBoundAtLevels(11, 4)

