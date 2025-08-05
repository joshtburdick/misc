#!/usr/bin/env python3
# Number of gates required for various things,
# for different sorts of gates.

import argparse
import itertools
import math
import pdb
import sys

import numpy as np

# note that comb() returns a float by default;
# for loop bounds, it needs the "exact=True" option,
# so that it returns an int
from scipy.special import comb

class UnboundedFanInNandBasis:
    """Bounds on number of unbounded-fan-in NAND gates."""

    def __init__(self):
        pass

    def num_functions(self, num_inputs, max_gates):
        """Number of functions implementable with some number of gates.

        This is an upper bound (for various reasons), but should still
            allow lower-bounding the number of gates.
        num_inputs: the number of inputs to the function
        max_gates: the maximum number of gates to include
        Returns: a numpy vector f with shape (max_gates+1,), such that
            f[i] is "number of functions using exactly i gates"
        """
        f = np.full(max_gates+1, None)
        # ??? this is the constant 1? or 0?
        f[0] = 1
        # note that this allows a gate with no inputs, which
        # presumably outputs 1
        for g in range(1, max_gates+1):
            f[g] = (2**(num_inputs+(g-1))) * f[g-1]
        # on second thought, assuming no circuits have zero gates
        f[0] = 0
        return f

    def expected_gates(self, num_inputs, lg_num_functions):
        """Lower bound on E[ # gates ] to compute some # functions.

        Note that this won't, in general, be an integer.
        num_inputs: the number of inputs to the circuit
        lg_num_functions: the number of functions (on log_2 scale);
            this can be an np.array
        """
        b = num_inputs - 0.5
        # the "-1" here is because this is the average, not the max.
        g = np.sqrt(2*lg_num_functions + b*b) - b - 1
        return g

    def and_upper_bound(self):
        """Upper bound on computing AND of two functions."""
        return 2

    def or_upper_bound(self):
        """Upper bound on computing OR of two functions.

        If A and B are the sets of cliques, then this is
        relative to the sum of the sizes of the two
        corresponding circuits.
        """
        # For unbounded fan-in, I think we save one gate (because we can
        # combine all the wires going into each output gate, feed them into
        # one of the output gates, and discard the other)
        return -1 

    def not_upper_bound(self):
        """Upper bound on computing NOT."""
        return 1

    def or_of_and_upper_bound(self, num_and, num_or):
        """Upper bound of ORing some ANDs together.

        This computes the AND of `num_and` inputs, `num_or` times, and
        then ORs the results together.
        """
        return num_or + 1

    def zonked_gates(self):
        """Number of gates zonked by feeding in a zero to one vertex."""
        # ??? is this needed?
        return 1.

class TwoInputNandBasis:
    """Bounds on number of two-input NAND gates.

    ??? include basis with two-input "any function" gates?
    """
    def __init__(self):
        """Constructor."""
        pass

    def num_functions(self, num_inputs, max_gates): 
        """Number of functions implementable with some number of gates.

        Again, this is an upper bound.
        I think that this is essentially the counting described in
            Aaronson, Scott. P != NP, pp. FIXME
        num_inputs: the number of inputs to the function
        max_gates: the maximum number of gates to include
        Returns: a numpy vector f with shape (max_gates+1,), such that
            f[i] is "number of functions using exactly i gates"
        """
        f = np.full(max_gates+1, None)
        # ??? how to define this?
        f[0] = 0
        # ??? this is the constant 0?
        f[1] = 1
        for g in range(2, max_gates+1):
            # each additional gate can use:
            # - an input
            # - a previous gate
            # - the constant 1 (thus changing a two-input NAND to a NOT)
            f[g] = comb(num_inputs + (g-1) + 1, 2, exact=True) * f[g-1]
        # FIXME double-check this
        return f

    def and_upper_bound(self):
        """Upper bound on computing AND of two functions."""
        # we can use a NAND, and then invert the output
        return 2
       
    def or_upper_bound(self):
        """Upper bound on computing OR of two functions.

        Returns: upper bound on number of gates to compute OR of two circuits,
            relative to their sum.
        """
        # invert the two inputs, then use a NAND
        return 3

    def not_upper_bound(self):
        """Upper bound on negating a function."""
        return 1

    def or_of_and_upper_bound(self, num_and, num_or):
        """Upper bound of ORing some ANDs together.

        This computes the AND of `num_and` inputs, `num_or` times, and
        then ORs the results together.
        """
        gates_for_and = 2 * (num_and-1)
        # ??? I think that, since the AND circuits all end with inverters,
        # we can just drop the final inverters, and use NAND gates to compute OR ?
        return num_or * gates_for_and + (num_or - 1)

    def zonked_gates(self):
        """Number of gates zonked by feeding in a zero to one vertex."""
        # ??? is this needed?
        return 1.

