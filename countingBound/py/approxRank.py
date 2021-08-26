#!/usr/bin/env python3
# Estimates rank of various things.

import matplotlib
import matplotlib.pyplot as plt

import numpy
import scipy.stats
import scipy.special

def rankBound(n, k):
    """Estimate of rank for finding some number of cliques.

    n: size of the input graph
    k: size of the cliques to find
    Returns: a vector r of length N+1, where N = choose(n, k),
        such that r_i is an estimate of the rank (in the list of all
        circuits) of the smallest circuit finding i cliques.
    """
    # number of cliques
    numCliques = scipy.special.comb(n, k, exact=True)
    # the bound (initially all 0)
    r = numpy.zeros([numCliques+1])
    # probability of a clique being "missed" when you feed in a 0
    # to a given edge
    p = 1 - (scipy.special.comb(k,2) / scipy.special.binom(n,2))
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
        r[i] = sum(w * r[0:i]) + scipy.special.comb(numCliques, i) / 2
    return r

b = rankBound(6,3)

plt.plot(range(21), b)
plt.title('n = 6, k = 3')
plt.xlabel('Level')
plt.ylabel('Rank')
# force x-axis to be plotted as integers
# ax = plt.figure().axes
# ax.xaxis.get_major_locator().set_params(integer=True)
plt.gca().xaxis.get_major_locator().set_params(integer=True)
plt.savefig('rank.png')

