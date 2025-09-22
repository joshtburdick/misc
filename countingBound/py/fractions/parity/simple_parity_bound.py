#!/usr/bin/env python3
# Further attempt at simplified bound for clique parity.
# FIXME
# - ideally `max_gates` should be computed automatically

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

class SimpleParity:
    """Attempt at a simpler bound for the clique parity problem.
    """

    def __init__(self, n, k, max_gates):
        """Constructor.

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        max_gates: maximum number of gates to include
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        self.max_gates = max_gates
        self.num_input_edges = comb(n, 2)
        # number of possible cliques
        self.num_possible_cliques = comb(n, k)
        # ??? make this an arg?
        self.basis = gate_basis.TwoInputNandBasis()
        # vector of counts of functions, using the naive upper bound
        self.function_count_low = FIXME
        # vector of counts of functions, using the naive upper bound,
        # _and_ a circuit for CLIQUE-parity for _all_ the cliques;
        # that circuit has index 0 in this vector
        self.function_count_high = FIXME





    def get_function_counts(self, upper_half):
        """Gets function counts per gate, for 0..N/2 cliques.

        Returns: a vector of length"""


    def get_function_counts_large_using_all(self):
        """Gets function counts, for N/2+1..N gates, using a CLIQUE-PARITY circuit.


        """
        FIXME


    def try_bound(self, max_gates):
        """Tests whether a circuit for clique parity would run into the counting bound.


        """
        # Get the counting bound for the number of possible functions constructable
        # using _up to_ a given number of gates.
        num_functions = self.basis.num_functions(comb(self.n, 2), self.max_gates+1).cumsum()



        pdb.set_trace()
        


def parse_args():
    parser = argparse.ArgumentParser(
        description='Bounds circuit size for detecting subsets of cliques.')
    parser.add_argument("n", type=int,
        help="Number of vertices")
    parser.add_argument("k", type=int,
        help="Size of clique")
    parser.add_argument("max_gates", type=int,
        help="Maximum number of gates to consider")
    parser.add_argument("--dump-lp",
        help="Dump LP problem statement to a file",
        action="store_true")
    parser.add_argument("--max-gates-search",
        help="Search for smallest number of gates which doesn't crash",
        action="store_true")
    parser.add_argument("--result-file",
        help="Write result to indicated file (rather than stdout)")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    b = SimpleParity(args.n, args.k, args.max_gates)


    for i in range(1, 5):
        print(b.try_bound(i))
        


#    if args.result_file:
#        with open(args.result_file, "wt") as f:
#            bounds.to_csv(f, index=False)

