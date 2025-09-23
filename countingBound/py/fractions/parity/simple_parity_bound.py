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
        self.function_counts_low = self.get_function_counts_low()
        # vector of counts of functions, using the naive upper bound,
        # _and_ a circuit for CLIQUE-parity for _all_ the cliques
        self.function_counts_high = self.get_function_counts_high_using_all()

    def get_function_counts_low(self):
        """Gets function counts per gate, for 0..N/2 cliques.

        """
        counts = np.zeros(self.max_gates+1, dtype=object)
        # 0 is the (hypothetical) number of gates in CLIQUE-PARITY
        counts[0] = comb(N, N)    # == 1; just emphasizing the reasoning here :)
        # i is the number of cliques included
        # i here will be the number of cliques being _removed_ from all of them,
        # by XORing with CLIQUE-PARITY.
        for i in range(1, self.N//2):
            # we need to check for i cliques, and XOR them together
            num_gates = i*self.gate_basis.and_cost(self.num_input_edges) + (i-1)*self.gate_basis.xor_cost()
            if num_gates > self.max_gates:
                return counts
            counts[ num_gates ] = comb(N, self.N//2 - i)

    def get_function_counts_high_using_all(self):
        """Gets function counts, for N/2+1..N cliques, using a CLIQUE-PARITY circuit.

        Returns: a vector of length max_gates+1, giving the number of functions
            constructable using the naive bound, and a CLIQUE-PARITY circuit.
            Here, index 0 corresponds to the size of the CLIQUE-PARITY circuit;
            so this will need to be offset, based on the conjectured size
            of the CLIQUE-PARITY circuit.
        """
        counts = np.zeros(self.max_gates+1, dtype=object)
        # the number of functions with as many gates in CLIQUE-PARITY (hypothetically)
        counts[0] = 1
        # i here will be the number of cliques being _removed_ from all of them,
        # by XORing with CLIQUE-PARITY.
        for i in range(1, self.N//2):
            # we need to check for i cliques, and XOR them together, _and_ with CLIQUE-PARITY
            num_gates = i*self.gate_basis.and_cost(self.num_input_edges) + i*self.gate_basis.xor_cost()
            if num_gates > self.max_gates:
                return counts
            counts[ num_gates ] = comb(N, i)
        return counts

    def try_bound(self, num_gates):
        """Tests whether a circuit for clique parity would run into the counting bound.

        num_gates_in_clique_parity: the assumed number of gates in a circuit for
            CLIQUE-PARITY
        Returns: the minimum amount by which the number of CLIQUE-PARITY functions
            "clears the ceiling" of the number of constructable functions (with a
            given number of gates.) If this is negative, then (for some number of gates),
            the number of CLIQUE-PARITY function exceeds the amount of "space"
            of functions constructable with that number of gates.
        """
        # Get the counting bound for the number of possible functions constructable
        # using _up to_ a given number of gates.
        num_constructable_functions = self.basis.num_functions(
            comb(self.n, 2), self.max_gates+1)
        num_constructable_functions_total = num_constructable_functions.cumsum()
        # Number of functions 
        num_functions = (self.function_counts_low()
            + np.pad(self.function_counts_high(), num_gates)[ : (self.max_gates+1) ])
        num_functions_total = num_functions.cumsum()

        # Find the minimum distance between these
        # ??? return the number of gates at which this is smallest?
        diff = num_constructable_functions_total - num_functions_total
        return np.min(diff)

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

