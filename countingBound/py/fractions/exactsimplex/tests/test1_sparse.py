#!/usr/bin/env python3
# test of sparse version, using problem from the original code

import fractions

# we assume this is being run from ../..
import sys
sys.path.append("./")

from exactsimplex.sparse import simplex

def list_to_fractions(x):
    return dict(enumerate(map(fractions.Fraction, x)))

def lol_to_fractions(x):
    return dict(enumerate(map(list_to_fractions, x)))


if __name__ == "__main__":
    c = list_to_fractions([300, 250, 450,  0,0,0,0])
    A = lol_to_fractions([
        [15, 20, 25,  1,0,0,0],
        [35, 60, 60,  0,1,0,0],
        [20, 30, 25,  0,0,1,0],
        [0, 250, 0,  0,0,0,-1]])
    b = list_to_fractions([1200, 3000, 1500, 500])

    # add slack variables by hand
    # A[0] += [1,0,0,0]
    # A[1] += [0,1,0,0]
    # A[2] += [0,0,1,0]
    # A[3] += [0,0,0,-1]  # is this the objective function?
    # c += [0,0,0,0]

    t, s, v = simplex(c, A, b)
    print(s)
    print(v)

