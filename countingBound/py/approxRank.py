#!/usr/bin/env python3
# Estimates rank of finding various numbers of functions.

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
    for i in range(k+1, n+1):
        # loop through the number of cliques with that many vertices
        maxCliques = comb(i, k, exact=True)
        for j in range(maxCliques+1):
            # the number of cliques which include an arbitrary vertex
            numCliquesIncludingVertex = comb(i-1, k-1, exact=True)
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
            r[(i,j)] = (
                    # the expected rank of functions which are definitely smaller
                    # (a weighted sum of the number of cliques "left over" after
                    # zonking cliques from one vertex) ...
                    sum( w * np.array([r[(i,j-z1)] for z1 in z]) )
                    # ... plus half the number of functions (with this many
                    # cliques, defined on this many vertices) being added
                    + comb(maxCliques, j) / 2)

    return r


##### plotting

n = 6
k = 3
maxCliques = comb(n, k, exact=True)

bound1 = rankBoundZeroedEdges(n, k)
bound2All = rankBoundZeroedVertices(n, k)
bound2 = [bound2All[n,k] for k in range(maxCliques+1)]

plt.plot(range(maxCliques+1), bound1, label="Edges")
plt.plot(range(maxCliques+1), bound2, label="Vertices")

plt.title('n = ' + str(n) + ', k = ' + str(k))
plt.xlabel('Level')
plt.ylabel('Rank')
# force x-axis to be plotted as integers
# ax = plt.figure().axes
# ax.xaxis.get_major_locator().set_params(integer=True)
plt.gca().xaxis.get_major_locator().set_params(integer=True)
plt.savefig('rank.png')

