#!/usr/bin/env python3
# test, from the original code

# we assume this is being run from ../..
import sys
sys.path.append("./")

from exactsimplex.simplex import simplex

if __name__ == "__main__":
   c = [300, 250, 450]
   A = [[15, 20, 25], [35, 60, 60], [20, 30, 25], [0, 250, 0]]
   b = [1200, 3000, 1500, 500]

   # add slack variables by hand
   A[0] += [1,0,0,0]
   A[1] += [0,1,0,0]
   A[2] += [0,0,1,0]
   A[3] += [0,0,0,-1]
   c += [0,0,0,0]

   t, s, v = simplex(c, A, b)
   print(s)
   print(v)

