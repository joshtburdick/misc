#!/usr/bin/env python3
# Estimates rank of various functions.

import pdb

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb, binom
from scipy.stats import hypergeom

def rankBound(n, k):
    """Estimate of rank for finding some number of cliques.

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

b = rankBound(6,3)
pdb.set_trace()

plt.plot(range(21), b)
plt.title('n = 6, k = 3')
plt.xlabel('Level')
plt.ylabel('Rank')
# force x-axis to be plotted as integers
# ax = plt.figure().axes
# ax.xaxis.get_major_locator().set_params(integer=True)
plt.gca().xaxis.get_major_locator().set_params(integer=True)
plt.savefig('rank.png')

