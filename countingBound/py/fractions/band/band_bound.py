#!/usr/bin/env python3
# Bound using a "band" of levels.

import argparse
import fractions
import itertools
import math
import pdb
import sys

import matplotlib as plt
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
    max_band_width: when bounding level i, for w in [0,max_band_width),
        group together levels i..i+w.
    Returns: an array of shape [max_band_width, comb(n, k)], of
        lower bounds.
    """
    basis = gate_basis.UnboundedFanInNandBasis()
    num_inputs = comb(n, 2)
    N = comb(n, k)
    b = np.zeros([max_band_width, N+1])
    num_functions = [comb(N, i) for i in range(N+1)]
    for j in range(N+1):
        for w in range(max_band_width):
            total_functions = sum(num_functions[j:j+w+1])
            b[w,j] = max(0, basis.expected_gates(num_inputs, math.log2(total_functions)) - w)
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
    print(np.max(b, axis=1))

