#!/usr/bin/python3
# Estimate of maximal clique size.

import math
from scipy.special import comb

# number of vertices
n = 4
# size of edges
r = 2

# note that n = 4, r = 2 is fairly near the border of what's convenient for
# "writing out all the cases"

for k in range(r, n + 1):
    # odds of a k-edge clique occurring
    probClique = 2 ** (-comb(k, r))
    # odds of it being covered by a given other vertex
    probOtherCover = 2 ** (-comb(k, r-1))
    # odds of it being "missed" by all the other vertices
    probAllMiss = (1 - probOtherCover) ** (n - k)
    # I think this should be the total number of maximal cliques of
    # some size, across all the possible hypergraphs
    # (or an estimate of that)
    numCliques = (2 ** comb(n,r)) * comb(n,k) * probClique * probAllMiss
    print(' '.join([str(n), str(r), str(k), str(numCliques)]))

# print(comb(10,4))

