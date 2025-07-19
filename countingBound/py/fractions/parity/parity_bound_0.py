#!/usr/bin/env python3
# Bound for clique parity, without using an LP.

import argparse
import fractions
import itertools
import math
import pdb
import sys

sys.path.append("..")   # XXX

import numpy as np
import pandas
import scipy.special
import scipy.stats

import gate_basis

# Wrapper for comb(), with exact arithmetic.
def comb(n, k):
    return scipy.special.comb(n, k, exact=True)

# Hypergeometric distribution, returning an exact fractions.Fraction .
def hyperg_frac(N, K, n, k):
    # based on https://en.wikipedia.org/wiki/Hypergeometric_distribution
    # note that we don't try to optimize this
    return fractions.Fraction(
        comb(K, k) * comb(N-K, n-k),
        comb(N, n))

class CliqueParity:
    """Attempt at bound for the clique parity problem.
    """

    def __init__(self, n, k, max_gates):
        """Constructor.

        ??? Should the "max number of gates considered" just be
        determined by the upper bound, based on n and k?)

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        max_gates: maximum number of gates to consider
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        self.max_gates = max_gates

        self.num_inputs = comb(n, 2)
        self.num_possible_cliques = comb(n, k)

        self.basis = gate_basis.TwoInputNandBasis()

        # Precompute lower bound on number of functions
        # (indexed by number of gates).
        self.num_possible_functions = self.basis.num_functions(
            self.num_inputs, self.max_gates)

        # Precompute upper bound on number of functions which compute
        # clique parity, up to the number of gates we're considering.
        self.upper_bound = np.zeros(self.max_gates+1, dtype=object)
        for i in range(self.num_possible_cliques // 2):
            num_gates = self.num_gates_upper_bound(i)
            if num_gates > self.max_gates:
                break
            self.upper_bound[ num_gates ] = comb(self.num_possible_cliques, i)

    def check_bound(self, num_gates_for_clique_parity):
        """Checks whether the given bound for clique parity is feasible.

        num_gates_for_clique_parity: the number of gates for clique parity
            which we're testing if it's feasible
        Returns: the minimum (across all number of gates) of the
            number of functions available, minus the number of
            functions needed. (If this is negative, then clique parity
            requires more gates.)
        """
        # lower bound on number of functions needed to compute
        # clique parity for "large" sets of cliques
        large_bound = np.concat([
            np.zeros(num_gates_for_clique_parity + 4, dtype=object),
            self.upper_bound])
        large_bound = large_bound[ : (self.max_gates+1) ]
        pdb.set_trace()
        diff = self.num_possible_functions - (self.upper_bound + large_bound)
        return np.min(diff)

    def num_gates_upper_bound(self, num_cliques):
        """Naive upper bound on computing parity of some number of cliques."""
        if num_cliques == 0:
            return 1
        # number of gates to detect each individual clique, by ANDing the edges
        gates_for_cliques = num_cliques * 2 * (comb(k,2)-1)
        # number of gates to XOR these together
        gates_for_xor = 4 * (num_cliques-1)
        return gates_for_cliques + gates_for_xor


def get_bound(n, k, max_gates):
    """Gets bound for some problem size.

    n, k, max_gates: problem size
    """
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, max_gates={max_gates}]\n')
    bound = CliqueParity(n, k, max_gates)
    for i in range(10):   # max_gates):
        slack = bound.check_bound(i)
        print(f"i={i}  {slack}")
#        if slack >= 0:
#            return i
    return 0

def parse_args():
    parser = argparse.ArgumentParser(
        description='Bounds circuit size for detecting subsets of cliques.')
    parser.add_argument("n", type=int,
        help="Number of vertices")
    parser.add_argument("k", type=int,
        help="Size of clique")
    parser.add_argument("max_gates", type=int,
        help="Maximum number of gates to consider")
    parser.add_argument("--result-file",
        help="Write result to indicated file (rather than stdout)")
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    n = args.n
    k = args.k
    max_gates = args.max_gates

    bound = get_bound(n, k, max_gates)
    print(f"bound={bound}")

    sys.exit(0)
    if args.result_file:
        with open(args.result_file, "wt") as f:
            bounds.to_csv(f, index=False)
    else:
        bounds.to_csv(sys.stdout, index=False)

