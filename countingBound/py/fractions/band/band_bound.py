#!/usr/bin/env python3
# Bound using a "band" of levels.

import argparse
import fractions
import itertools
import math
import pdb
import sys

import numpy as np
import pandas
import scipy.special
import scipy.stats

sys.path.append("..")   # XXX
        
import gate_basis

# Wrapper for comb(), with exact arithmetic.
def comb(n, k):
    return scipy.special.comb(n, k, exact=True)

def band_bound(n, k, max_band_width):
    """Computes the 'band bound'.

    n, k: Problem size.
    max_band_width: when bounding level i, for j in [0,max_band_width),
        group together levels i..i+j.
    Returns: an array of shape [max_band_width, comb(n, k)], of
        lower bounds.
    """
    basis = gate_basis.UnboundedFanInNandBasis()
    num_inputs = comb(n, 2)
    N = comb(n, k)
    b = np.zeros([max_band_width, N], dtype=int)
    num_functions = comb(N, range(N+1))
    for i in range(N+1):
        for j in range(max_band_width):
            total_functions = sum(num_functions[i:i+j])
            b[i,j] = max(0, basis.expected_gates(num_inputs, log2(total_functions)) - j)
    return b

def parse_args():
    parser = argparse.ArgumentParser(
        description='Bounds circuit size, by grouping nearby numbers of cliques.'
            ' Writes results in currint working directory.')
    parser.add_argument("n", type=int,
        help="Number of vertices")
    parser.add_argument("k", type=int,
        help="Size of clique")
    parser.add_argument("max_band_width", type=int,
        help="Maximum 'band width' to consider")
   return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    b = band_bound(args.n, args.k, args.max_band_width)





