#!/usr/bin/env python3
# Linear programming bound on rank of different functions.
# This includes several constraints; not all of them may be useful.
# In this version, we zero out a vertex at a time.

import numpy as np
import scipy.optimize
from scipy.special import comb

class LpBound:
    """Computes a bound on the rank of finding cliques.

    Here, we track the number of cliques, and the total number of
    vertices used by all of the cliques (which hopefully will be
    easier to deal with.)
    """

    def __init__(self, n, k):
    """ Constructor
    n: number of vertices
    k: clique size
    """
    self.n = n
    self.k = k
    # mapping from tuples (numVertices, numCliques) to index in the LP
    self.varIndex = {}
    index = 0
    # set up the mapping
    for i in range(k, n+1):
        for j in range(0, comb(i, k)):
            self.varIndex[(i,j)] = index
            index += 1
    # this is the total number of variables we're solving for
    self.numVariables = index - 1
    # these store the constraints, as lists (for easy appending,
    # since it's not clear how many there will be).
    # A is stored as a list of triples
    # (numVertices, numCliques, coefficient), which hopefully
    # will be easier to understand...
    self.A = []
    # ... and b as a list of numbers
    self.b = []
    # the bounds matrices (initially undefined)
    self.A_ub = None
    self.b_ub = None



    def addConstraint(self, A, b):
        """Adds one row to the constraints.

        A: a list of (numVertices, numCliques, coefficient) triples
        b: the corresponding lower bound
        Side effects: adds a row to the bound, of the form "Ax >= b"
        """
        self.A.append(A)
        self.b.append(b)

    def addVertexZeroConstraint(self):
        """Adds constraint from restricting some vertex's edges to 0.

        Note that punctuating the possessive of a word ending in 'x'
        is just problematic.
        Side effects: adds a constraint on expected.
        """
        # j is the number of vertices _after_ a vertex is zeroed out
        for j in range(self.k, self.n):
            # i is the number of cliques
            a = [(i,j,FIXME)] for i in range(   , )]


        pass

    def addVertexTotalConstraint(self):
        """Constraint on rank of sets with some number of vertices.

        Side effects: for each possible number of vertices, adds
        a constraint on the total rank of the sets with that
        many vertices.
        """
        # j is the number of vertices
        for j in range(self.k, self.n):
            # constraint on the "weighted average" of these



    def setBounds():
        """Sets the bounds matrices, A_ub and b_ub."""
        # if these are already computed, skip this
        if self.A_ub and self.B_ub:
            return
        # converts from "list of numbers" to a row of A
        def constraintRow(A_list):
            row = np.zeros(self.numVariables)
            for entry in A_list:
                (numVertices, numCliques, a) = entry
                row[ self.varIndex[(numVertices, numCliques)] ] = a
            return row 
        # convert A and b to np.array objects (note that both are
        # negated, since the solver is solving Ax <= b).
        self.A_ub = - np.stack([constraintRow(a1) for a1 in self.A])
        # b is just converted into a column vector
        self.b_ub = - np.array(self.b)
 
    def solve(numVertices):
        """Solves the linear system.
        
        Note that by default, the solver constrains all x >= 0,
        so we don't add that constraint.
        numVertices: the number of vertices to solve for
        Returns: a numpy array, of the minimum rank of finding
            all the cliques.
        """
        # the objective function: how low can finding all the
        # cliques go?
        c = np.zeros(self.numVariables)
        c[ self.numHittingEdge(), self.numNotHittingEdge() ] = 1
        



