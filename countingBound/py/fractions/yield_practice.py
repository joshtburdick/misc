#!/usr/bin/env python3
# Practicing using "yield".

import pdb

from scipy.special import comb

def function_counts(num_inputs, max_fan_in=2):
    """Gets counts of functions

    num_inputs: number of inputs to the circuit
    max_fan_in: maximum number of inputs (None for "unbounded")
    """
    num_gates = 0
    # we assume that with zero gates, we can either
    # "pass through" an input, or implement constants 0 and 1
    num_functions = num_inputs + 2
    # number of inputs to the next gate
    num_gate_inputs = num_inputs + 2
    while True:
        yield (num_gates, num_functions)
        # "pool" of inputs to new gates are the inputs,
        # outputs of previous gates, or 0 or 1
        num_gate_inputs = num_inputs + num_gates + 2
        num_gates += 1
        if max_fan_in:
            num_functions *= comb(num_gate_inputs,
                max_fan_in, exact=True)
        else:
            num_functions *= 2 ** num_gate_inputs

y = function_counts(6, max_fan_in=None)
print(next(y))
print(next(y))
print(next(y))

# pdb.set_trace()

