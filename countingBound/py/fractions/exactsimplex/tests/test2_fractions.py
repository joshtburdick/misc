#!/usr/bin/env python3
# another small test

import fractions

# we assume this is being run from ../..
import sys
sys.path.append("./")

from exactsimplex.simplex import simplex

def list_to_fractions(x):
    return list(map(fractions.Fraction, x))

def lol_to_fractions(x):
    return list(map(list_to_fractions, x))


if __name__ == "__main__":
   # minimizing sum of these
   c = list_to_fractions([-1, -1])
   A = lol_to_fractions([[1, 5], [2, 1]])
   b = list_to_fractions([-100, -100])

   # add slack variables by hand
   A[0] += [1,0]
   A[1] += [0,1]
   c += [0,0]

   t, s, v = simplex(c, A, b)
   print(t)
   print(s)
   print(v)

