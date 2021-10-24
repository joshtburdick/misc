#!/usr/bin/env python3
# Linear programming bound on rank of different functions.
# This includes several constraints; not all of them may be useful.
# In this version, we zero out a vertex at a time.
# Note that the nomenclature is confusing here.

import numpy as np
import pdb
import scipy.optimize
# note that comb() returns a float by default;
# for loop bounds, it needs the "exact=True" option,
# so that it returns an int
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
        # problem size
        self.n = n
        self.k = k
        # mapping from tuples (numVertices, numCliques) to
        # variable index in the LP
        self.varIndex = {}
        index = 0
        # set up the mapping of variable indices
        for i in range(k, n+1):
            for j in range(0, comb(i, k, exact=True)+1):
                self.varIndex[(i,j)] = index
                index += 1
        # this is the total number of variables we're solving for
        self.numVariables = index
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

    def addVertexTotalConstraint(self):
        """Constraint on rank of sets with some number of vertices.

        Note that the sets with more vertices also includes functions
        with a smaller number of vertices.

        Side effects: for each possible number of vertices, adds
        a constraint on the total rank of the sets with that
        many vertices.
        """
        # i is the number of vertices
        for i in range(self.k, self.n+1):
            # the number of possible cliques with that many vertices
            numCliques = comb(i, self.k, exact=True)
            # the number of functions with up to that many cliques
            numFunctions = 2 ** numCliques
            # constraint on the "weighted average" of these
            # (here, i is the number of cliques in the function)
            a = [(i, j, comb(numCliques, j) / numFunctions)
                    for j in range(0, numCliques+1)]
            # the weighted average should be at least
            # half the number of functions
            self.addConstraint(a, numFunctions / 2)

    def addVertexZeroExpectedConstraint(self):
        """Adds constraint from restricting some vertex's edges to 0.

        This constraint says that if you take a random graph with
        j+1 vertices, and zero out all the edges from one vertex,
        the rank of the resulting graph will be smaller.
        
        Note that punctuating the possessive of a word ending in 'x'
        is just problematic.
        ??? also add constraint that "zeroing out a vertex's edges
            strictly reduces rank"?
        Side effects: adds a constraint on expected rank.
        """
        A = []
        # i is the number of vertices _after_ a vertex is zeroed out
        # (and thus ranges up to n-1)
        for i in range(self.k, self.n):
            # maximum number of cliques which might be made
            # impossible, by zeroing out the edges connected to a vertex
            numCliquesZeroed = comb(i-1, self.k-1)
            # corresponding number of functions
            numFunctionsZeroed = 2 ** numCliquesZeroed
            # j is the number of cliques _after_ a vertex is zeroed out
            for j in range(0, comb(i, self.k)+1):
                # the rank, after a vertex is zeroed out
                a = [(i, j, 1.0)]
                # k is the number of cliques which were zeroed out
                # (this shouldn't throw a KeyError)
                a += [(i+1, j+k, -comb(maxCliquesZeroed, k))
                        for k in range(0, maxCliquesZeroed+1)]
                # the constraint is that "the expected rank after
                # zeroing out a clique is some amount higher than
                # the rank of what remains"
                b = 0
                self.addConstraint(a, b)

    def setBounds(self):
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
 
    def solve(self, numVertices):
        """Solves the linear system.
        
        Note that by default, the solver constrains all x >= 0,
        so we don't add that constraint.
        numVertices: the number of vertices in the function to minimize.
        Returns: a numpy array, of the minimum rank of finding
            all the k-vertex cliques in an m-vertex graph,
            for 0 <= m <= n. (Note that this is only minimized
            for m==numVertices; for the others, I'm curious what bounds
            it gives, but there's no guarantee that it's minimal.)
            For m < k, this is 0, but those cases are included
            for indexing convenience.
        """
        # set A_ub and b_ub (if they haven't been set already)
        self.setBounds()
        # the objective function: how low can the rank of finding
        # all the cliques (with that many vertices) be?
        c = np.zeros(self.numVariables)
        numCliques = comb(numVertices, self.k)
        c[ self.varIndex[(numVertices, numCliques)] ] = 1
        # solve
        r = scipy.optimize.linprog(c, self.A_ub, self.b_ub)
        return(r)

if __name__ == '__main__':
    print('in main')
    lp = LpBound(5,3)
    # this probably won't do much
    lp.addVertexTotalConstraint()
    # pdb.set_trace()
    r = lp.solve(5)
    print(r)
