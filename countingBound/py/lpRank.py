#!/usr/bin/env python3
# Linear programming bound on rank of different functions.
# This could almost certainly be simplified, but here we're
# trying to emphasize clarity.
# This includes several constraints; not all of them may be useful.

import numpy as np

import scipy.optimize
import scipy.special

class LpBound:
    """Computes a bound on the rank of finding cliques.

    n: number of vertices
    k: clique size
    """
    self.n = n
    self.k = k
    # the total number of functions
    self.numFunctions = scipy.special.comb(n, k)
    # the number of functions hitting an arbitrary edge,
    # (say, e_{12})
    self.numHittingEdge = scipy.special.comb(n-2, k-2)
    # number of variables in a constraint
    self.numVarsInConstraint = (2 * numFunctions
            + self.numHittingEdge *
            (self.numFunctions - self.numHittingEdge)
    # these store the constraints, as lists (for easy appending,
    # since it's not clear how many there will be).
    # ??? rename these?
    # A is stored as a list of numpy vectors...
    self.A = []
    # ... and b as a list of numbers
    self.b = []

    def getConstraintRow():
        """Gets a row of the constraint matrix A (of the right size)."""
        return np.array([ self.numVarsInConstraint ])

    def addConstraint(A, b):
        """Adds one row to the constraints.

        A: a numpy.array of constraint coefficients
        b: the corresponding bound
        Side effects: adds a row to the bound, of the form "Ax >= b"
        """
        # note that we negate both of these, since scipy.linsolve
        # inequalities are "Ax <= b"
        self.A.append(-A)
        self.b.append(-b)

    def C(self, i):
        """Index of C_i, E[rank(finding i cliques)]."""
        return i - 1

    def B(self, i):
        """Index of B_i,
        E[rank(finding i cliques which miss e_{12})].
        """
        return self.numFunctions + i - 1

    def A(self, i, j):
        """Index of A_{ij}, E[rank("higher-up functions")].

        Specifically, this is E[rank(finding a set of cliques,
        of which i miss e_{12}, and j hit e_{12}.
        """
        return (2 * self.numFunctions
                + i * self.numHittingEdge
                + j)

    def addBTotalConstraint():
        """Adds constraint on expected value of E[B].
       
        (It's possible it will make more sense to constrain E[C].)
        Side effects: adds the constraint.
        """
        a = np.zeros([ self.numVarsInConstraint ])
        # maximum number of cliques in B (if there are more,
        # some cliques must be in A)
        m = self.numFunctions - self.numHittingEdge
        # constrain the weighted average of all of these
        for i in range(m+1):
            a[ self.B(i) ] = scipy.special.comb(m, i) / (2 ** m)
        # what the weighted average should be
        self.addConstraint(a, 0.5 * (2 ** m))


    def addAConstraint():
        """Adds a constraint on A.

        This constraint 

        """
        pass

    def solve():
        """Solves the linear system.
        
        Note that by default, the solver constrains all x >= 0.

        Returns: a numpy array, of the minimum rank of C_i.
        """
        # convert A and b of the constraints into numpy objects
        A = np.stack(self.A, axis=FIXME)
        b = np.array(self.b)

        # the objective function: how low can C_N go?
        c = np.zeros([ self.numVarsInConstraint ])
        c[ self.C(self.numFunctions) ] = 1




