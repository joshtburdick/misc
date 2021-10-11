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
    # the total number of cliques
    self.numCliques = scipy.special.comb(n, k)
    # the number of sets of cliques (and functions) which
    # include an arbitrary edge, (say, e_{12})
    self.numHittingEdge = scipy.special.comb(n-2, k-2)
    # the number of sets of cliques which don't include e_{12}
    self.numNotHittingEdge = self.numCliques - self.numHittingEdge
    # these store the constraints, as lists (for easy appending,
    # since it's not clear how many there will be).
    # ??? rename these?
    # A is stored as a list of numpy matrices
    self.A = []
    # ... and b as a list of numbers
    self.b = []

    def getConstraintMatrix():
        """Gets a matrix corresponding to one constraint.
        
        The row index is the number of cliques in A (which
        intersect an arbitrarily-chosen edge), while
        the column index is the number of cliques in B.
        """
        return np.zeros([self.numHittingEdge+1,
            self.numNotHittingEdge+1])

    def addConstraint(A, b):
        """Adds one row to the constraints.

        A: a numpy.array of constraint coefficients
        b: the corresponding bound
        Side effects: adds a row to the bound, of the form "Ax >= b"
        """
        self.A.append(A)
        self.b.append(b)

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


    def addZeroRestrictionConstraint():
        """Adds effect of restricting e_{12} to 0.

        


        """
        pass

    def addPermutationConstraint():
        """Add constraint that "permuting vertices doesn't matter".

        XXX FIXME it's not clear what this should look like
        """
        pass

    def solve():
        """Solves the linear system.
        
        Note that by default, the solver constrains all x >= 0.

        Returns: a numpy array, of the minimum rank of each set
            of functions.
        """
        # convert A and b to np.array objects (note that both are
        # negated, since the solver is solving Ax <= b).
        A = - np.stack([a.reshape(FIXME) for a in self.A], axis=?)
        # b is just converted into a column vector
        b = - np.array(self.b)
        # the objective function: how low can finding all the
        # cliques go?
        c = self.getConstraintMatrix()
        c[ self.numHittingEdge(), self.numNotHittingEdge() ] = 1
        c = c.reshape(FIXME)




