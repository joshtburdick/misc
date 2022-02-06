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
    Note that this seems inaccurate, because "zonking" edges in this
    way means that we aren't sampling uniformly from circuits which
    find smaller numbers of cliques.
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

def zeroingStep(rankBoundA, sizeOfB):
    """Calculates a new bound, after zeroing out some edges.

    This is the "step case": given a bound for some set of cliques A,
        it gives a bound if some other set of cliques B were added.
    rankBoundA: numpy array (whose length implicitly gives the size
        of A). Its i'th entry is a lower bound on the expected rank
        of a circuit (picked uniformly randomly) which finds
        exactly i cliques.
    sizeOfB: number of cliques in B, which are being added.
        There should be an edge (or edges) which can be set to zero,
        which remove all the cliques in B, but leave all the cliques
        in A.
    Returns: a numpy array (like rankBoundA) whose i'th entry is
        a lower bound on the expected rank of a circuit which
        finds exactly i cliques, except in C = A union B.
    """
    sizeOfA = rankBoundA.shape[0]-1
    sizeOfC = sizeOfA + sizeOfB
    rankBoundC = np.zeros(sizeOfC+1)
    # loop through the possible combined sizes
    for numCliquesInC in range(sizeOfC+1):
        # bounds on how many of these cliques might be in A
        minCliquesInA = max(0, numCliquesInC - sizeOfB)
        maxCliquesInA = min(sizeOfA, numCliquesInC)
        # the possibilities for the number of cliques in A
        a = np.arange(minCliquesInA, maxCliquesInA+1)
        # the probability of some number of these being in A
        w = hypergeom(
            # number of possible cliques
            sizeOfC,
            # number of those present (in C = A union B)
            numCliquesInC,
            # number of possible cliques in A
            sizeOfA
            ).pmf(a)
        # the number of "new" functions we get by adding B
        # (which is the number in A union C, minus those
        # already counted in A)
        numNewFunctions = (comb(sizeOfC, numCliquesInC, exact=True)
                - comb(sizeOfA, numCliquesInC, exact=True))
        # combining these: the bound is a weighted sum of the
        # previous bound (from A), and half the number of "new" functions
        rankBoundC[ numCliquesInC ] = (w.dot( rankBoundA[ a ]) + numNewFunctions / 2)
    # return the bound, for A union B
    return rankBoundC

def rankBoundZeroedVertices1(maxVertices, k):
    """Estimate of rank for finding some number of cliques.

    This version feeds in zeros to all edges incident to a vertex.
    maxVertices: largest size of the input graph
    k: size of the cliques to find
    Returns: a dict r, keyed by numVertices,
        such that r[i] is a lower bound on the expected rank of
        the smallest circuit finding i cliques in a graph with
        <= numVertices vertices.
    """
    # initially there are two functions
    bound = { k: np.array([0, 1]) }
    # loop through the number of vertices
    for numVertices in range(k+1, maxVertices+1):
        bound[numVertices] = zeroingStep(bound[numVertices-1],
            comb(numVertices-1, k-1, exact=True))
    return bound

def rankBoundZeroingVertexEdges1(maxNumVertices, k, includePartialVertices=False):
    """Bounds function rank, zeroing a vertex' edges, one at a time.

    This version adds vertices, one at a time. For each, it adds
    edges one at a time, and tracks how many cliques would be lost
    by zeroing out that edge.
    XXX this function is getting looooong

    maxNumVertices: size of the input graph
    k: size of the cliques to find
    includePartialVertices: if True, include cases in which a vertex is only
        partially added. (If False, the (v,e) keys below will all have e==0.)
    Returns: a dict r, keyed by (v, e), where:
        v: the number of vertices "completely" added so far
            (which ranges up to maxNumVertices)
        e: the number of edges added (thus far) to the "new" vertex
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
        for numNewEdges in range(k-1, v):
            # The number of cliques created by adding an edge e;
            # all of these would be "zonked" if e were set to 0.
            # e's ends determine two vertices of these cliques, but
            # the remaining edges could be any subset of the remaining edges.
            numNewCliques = comb(numNewEdges-1, k-2, exact=True)
            # update the bound on rank for every possible number of cliques
            # (with the new edge, more cliques are possible)
            bound1 = zeroingStep(bound1, numNewCliques)
            # possibly also save the bound when the new vertex is
            # only partially added (mostly for debugging) ?
            if includePartialVertices and numNewEdges < v-1:
                bound[(v, numNewEdges)] = bound1
        # save the bound, for "no additional edges"
        # ??? shouldn't this be v+1 ?
        bound[(v,0)] = bound1
    return bound

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

def rankBoundZeroingVertexEdges(maxNumVertices, k, includePartialVertices=False):
    """Bounds function rank, zeroing a vertex' edges, one at a time.

    This version adds vertices, one at a time. For each, it adds
    edges one at a time, and tracks how many cliques would be lost
    by zeroing out that edge.
    XXX this function is getting looooong

    maxNumVertices: size of the input graph
    k: size of the cliques to find
    includePartialVertices: if True, include cases in which a vertex is only
        partially added. (If False, the (v,e) keys below will all have e==0.)
    Returns: a dict r, keyed by (v, e), where:
        v: the number of vertices "completely" added so far
            (which ranges up to maxNumVertices)
        e: the number of edges added (thus far) to the "new" vertex
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
        for numNewEdges in range(k-1, v):
            # number of possible cliques so far (formed from vertices
            # other than v, plus edges connected from v so far)
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
            # possibly also save the bound when the new vertex is
            # only partially added (mostly for debugging) ?
            if includePartialVertices and numNewEdges < v-1:
                bound[(v, numNewEdges)] = bound2
        # save the bound, for "no additional edges"
        # ??? shouldn't this be v+1 ?
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
    # compute various bounds, and plot them
    plt.figure()
    plt.plot(range(maxCliques+1),
        [comb(maxCliques, c)/2 for c in range(maxCliques+1)],
        label='Naive counting bound', alpha=0.5, linewidth=5)

    bound1 = rankBoundZeroedEdges(n, k)
    plt.plot(range(maxCliques+1), bound1, label='Zeroing out edges', alpha=0.6)

    # bound2All = rankBoundZeroedVertices(n, k)
    # bound2 = [bound2All[n,k] for k in range(maxCliques+1)]
    # plt.plot(range(maxCliques+1), bound2[n], label='Zonking vertices (take 2)', alpha=0.6)

    bound2a = rankBoundZeroedVertices1(n, k)
    plt.plot(range(maxCliques+1), bound2a[n], label='Zonking vertices (take 2)', alpha=0.6)

    bound3 = rankBoundZeroingVertexEdges1(n+1, k)
    bound3 = bound3[(n, 0)]
    plt.plot(range(maxCliques+1), bound3,
            label='Zonking vertices, an edge at a time')

    # bound of "half of all functions of size <= this"
    # FIXME rename?
    # also currently omitted, as it looks weirdly high
    if False:
        plt.plot(range(maxCliques+1),
            np.cumsum(np.array(
                [comb(maxCliques, c)/2 for c in range(maxCliques+1)])),
            label='"Half of all the functions"')
    plt.title('n = ' + str(n) + ', k = ' + str(k))
    plt.xlabel('Number of cliques')
    plt.ylabel('E[rank of functions]')
    plt.legend()
    # force x-axis to be plotted as integers
    plt.gca().xaxis.get_major_locator().set_params(integer=True)
    # plt.savefig('rank_n=' + str(n) + '_k=' + str(k) + '.pdf')
    plt.savefig('rank_n=' + str(n) + '_k=' + str(k) + '.png')

plotBoundAtLevels(6, 3)
if True:
    for n in range(6, 12):
        plotBoundAtLevels(n, 3)
        plotBoundAtLevels(n, 4)
        plotBoundAtLevels(n, 5)
if False:
    for k in range(6, 7):
        plotBoundAtLevels(2*k, k)

# this runs, but gives a not-great bound
if False:
    z = rankBoundZeroingVertexEdges(8, 3, includePartialVertices=True)
    for i in z.items():
        print(i[0], i[1].shape[0]-1)

if False:
    plotBoundAtLevels(6, 3)
    plotBoundAtLevels(11, 4)

