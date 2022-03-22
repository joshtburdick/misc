#!/usr/bin/env python3
# Linear programming bound on rank of different functions.
# This includes several constraints; not all of them may be useful.

import numpy as np
import pdb
import scipy.optimize
# note that comb() returns a float by default;
# for loop bounds, it needs the "exact=True" option,
# so that it returns an int
from scipy.special import comb
from scipy.stats import hypergeom

class LpBound:
    """Computes a bound on the rank of finding cliques.

    This version tries to do bookkeeping by zeroing out a distinguished edge, e.
    ??? should I rename the edge which is zeroed out?
    """

    def __init__(self, n, k):
        """ Constructor.

        ??? should I rename these?
        n: number of vertices
        k: size of clique we're looking for
        """
        # problem size
        self.n = n
        self.k = k
        # number of cliques possible
        self.max_cliques = comb(n, k, exact=True)
        # number of cliques which could be zeroed out when edge e is zeroed out
        self.max_cliques_zeroed = comb(n-2, k-2, exact=True)
        # how many cliques could be left over
        self.max_cliques_remaining = self.max_cliques - self.max_cliques_zeroed
        # mapping from tuples (numVertices, numCliques) to
        # variable index in the LP
        self.var_index = {}
        # set up the mapping of variable indices
        for i in range(self.max_cliques_zeroed+1):
            for j in range(self.max_cliques_remaining+1):
                self.var_index[(i,j)] = len(self.var_index)
        # these store the constraints, as lists (for easy appending,
        # since it's not clear how many there will be).
        # A is stored as a list of numpy arrays
        self.A = []
        # ... and b as a list of numbers
        self.b = []

    def add_constraint(self, A, b):
        """Adds one row to the constraints.

        A: a list of (index, coefficient) pairs, where "index" is
            a key in varIndex
        b: the corresponding lower bound
        Side effects: adds a row to the bound, of the form "Ax >= b"
        """
        # converts from "list of numbers" to a row of A
        A_row = np.zeros(len(self.var_index))
        for entry in A:
            (i, a) = entry
            A_row[ self.var_index[ i ] ] = a
        self.A.append(A_row)
        self.b.append(b)

    def add_counting_bound_constraints(self):
        """Adds constraints based on the counting bound.

        For each "level" of "total number of cliques found", this
        adds a bound, based on the counting bound.
        """
        # loop through the number of cliques
        for num_cliques in range(self.max_cliques+1):
            # the maximum number of cliques containing edge e
            # (these won't actually be zeroed)
            max_cliques_zeroed = min(num_cliques, self.max_cliques_zeroed)
            # the probability of some number of cliques containing edge e
            h = hypergeom(
                # number of possible cliques
                self.max_cliques,
                # number of those present
                num_cliques,
                # number of cliques which could intersect edge e
                max_cliques_zeroed)
            # FIXME explain this better
            # first, the case when no cliques intersect edge e
            A = [((0, num_cliques), h.pmf(0) - 1)]
            # here, z is the number of cliques which _do_ intersect edge e
            A += [((z, num_cliques-z), h.pmf(z))
                for z in range(1, max_cliques_zeroed+1)]
            # the bound is half the number of functions
            b = comb(self.max_cliques, num_cliques, exact=True) / 2
            self.add_constraint(A, b)

    def add_edge_zeroing_constraints(self):
        """Adds constraints based on zeroing out an edge.

        The bound is that all the sets of cliques which intersect
        edge e have a higher expected rank than the sets
        remaining after feeding in a 0 to edge e.
        """
        # number of functions with a clique overlapping edge e, plus
        # 1 (for the function without any cliques overlapping e)
        num_higher_functions = 2 ** self.max_cliques_zeroed
        # the probabilities of some number of cliques being zeroed
        p = comb(self.max_cliques_zeroed, range(self.max_cliques_zeroed+1))
        # loop through the number of cliques "left over"
        pdb.set_trace()
        for j in range(self.max_cliques_remaining):
            # this has a "-1" because the functions in which some cliques
            # include edge e are all higher than the functions in which
            # no cliques include edge e
            A = [(0, p[0] - 1)]
            # constraints in the case of functions including a clique
            # which includes edge e
            A += [((i,j), p[i])
                for i in range(1, self.max_cliques_zeroed+1)]
            b = num_higher_functions / 2
        self.add_constraint(A, b)

    def solve(self):
        """Solves the linear system.
        
        Note that by default, the solver constrains all x >= 0,
        so we don't add that constraint.
        FIXME: add option to minimize finding only some number
            of cliques (rather than all of them) ?
        Returns: the LP result
        """
        # convert A and b to np.array objects (note that both are
        # negated, since the solver is solving Ax <= b).
        self.A_ub = - np.stack(self.A)
        # b is just converted into a column vector
        self.b_ub = - np.array(self.b)
        # the objective function: how low can the rank of finding
        # all the cliques (with that many vertices) be?
        c = np.zeros(self.numVariables)
        c[ self.varIndex[(self.max_cliques_zeroed, self.max_cliques_remaining)] ] = 1
        # solve
        r = scipy.optimize.linprog(c, self.A_ub, self.b_ub)
        # for now, we return the entire result (rather than just
        # the result), in case it's useful for debugging
        return r

if __name__ == '__main__':
    print('in main')
    lp = LpBound(5,3)
    pdb.set_trace()
    lp.add_counting_bound_constraints()
    lp.add_edge_zeroing_constraints()
    r = lp.solve()
    print(r)

