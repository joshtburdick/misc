#!/usr/bin/env python3
# Attempt based on zeroing out vertices.
# This is derived from ../lp_gate_bound_6.py, which
# was using the LP relaxation of an IP.
# This is intended to allow using an IP, or the LP relaxation.

import argparse
import fractions
import itertools
import math
import pdb
import sys

import numpy as np
import pandas
import scipy.special
import scipy.stats

import gate_basis
import scip_helper

# Wrapper for comb(), with exact arithmetic.
def comb(n, k):
    return scipy.special.comb(n, k, exact=True)

# Hypergeometric distribution, returning an exact fractions.fraction .
def hyperg_frac(N, K, n, k):
    # based on https://en.wikipedia.org/wiki/Hypergeometric_distribution
    # note that we don't try to optimize this
    return fractions.Fraction(
        comb(K, k) * comb(N-K, n-k),
        comb(N, n))

class LpVertexZeroing:
    """Attempt at bound by zeroing out vertices.
    """

    def __init__(self, n, k, max_gates):
        """Constructor gets graph info, and sets up variable names.

        This will have two groups of variables, labeled either
            by 2-tuples or 3-tuples.
        The 2-tuples will be of the form ("x", v, c) where:
            - "v" is the number of vertices, with k <= v <= n
            - "c" is the number of cliques, with 0 <= c <= {n choose v}
            - Each variable will be the expected number of gates in
            the sets of cliques with _up to_ v vertices.
        The 3-tuples will have labels ("g", v, c, gates) where:
            - "v" is as above
            - "c" is as above
            - "gates" is the number of gates
            - Each variable will be the number of sets of cliques
              (with size given by v and c), with smallest circuit
              having "gates" gates.
        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        max_gates: maximum number of gates to include
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        self.max_gates = max_gates

        self.expected_num_gates_vars = []
        self.num_gates_dist_vars = []
        for v in range(k, n+1):
            for c in range(0, comb(v, k) + 1):
                # variables for expected number of gates
                self.expected_num_gates_vars += [(v, c)]
                for g in range(max_gates+1):
                    # variables for counts of numbers of functions
                    # with some number of gates
                    self.num_gates_dist_vars += [(v, c, g)]

        # wrapper for LP solver
        self.lp = scip_helper.SCIP_Helper(
            self.expected_num_gates_vars + self.num_gates_dist_vars)

        # self.basis = gate_basis.UnboundedFanInNandBasis()
        self.basis = gate_basis.TwoInputNandBasis()

    def add_level_constraints(self):
        """Adds constraints on functions at some "level".

        By "level", we mean "number of vertices" _and_ "number of cliques".

        The constraints (both of which are equalities) are:
        - on the total number of functions at that "level", and
        - connecting the counts with that "level"'s expected gate count
        """
        # loop through number of vertices
        for v in range(self.k, self.n+1):
            num_possible_cliques = comb(v, self.k)
            for c in range(num_possible_cliques+1):
                num_functions = comb(num_possible_cliques, c)
                # add constraint that these sum to the number of functions
                self.lp.add_constraint(
                    [((v, c, i), 1.) for i in range(self.max_gates+1)],
                    '=', num_functions)
                # add constraint defining expected number of gates
                A = [((v, c, i), i)
                    for i in range(1, self.max_gates+1)]
                self.lp.add_constraint(A + [((v, c), -num_functions)],
                    '=', 0.)

    def add_counting_bound(self):
        """Adds counting bounds, for a given number of gates.

        This is a lower bound on the number of functions with
        some number of gates.
        """
        # we have a separate bound for each possible number of vertices
        # (which will affect the number of inputs, and thus the number
        # of possible functions)
        for v in range(self.k, self.n+1):
            # number of possible cliques
            num_possible_cliques = comb(v, self.k)
            # number of possible functions for each possible number of gates
            # (with number of inputs based on number of vertices)
            num_functions = self.basis.num_functions(comb(v, 2), self.max_gates)
            # upper-bound "total number of functions with this many gates"
            for g in range(self.max_gates+1):
                self.lp.add_constraint(
                    [((v, c, g), 1.)
                        for c in range(num_possible_cliques+1)],
                    '<=', num_functions[g])

    def add_vertex_zeroing_constraints(self):
        """Adds constraints from zeroing out vertices.

        Note that this only interacts with the expected gate count
        variables (that is, the 2-tuples), and so is the same as before.

        The sets of cliques are:
        C: the cliques in the larger set, using (up to) k+1 vertices
        B: the cliques zeroed by feeding in zeros to a vertex of C
        A: the cliques which are left over
        """
        # loop through number of vertices in larger graph
        for v in range(self.k+1, self.n+1):
            # loop through number of cliques in that graph
            # (note that the trivial bound when C_size==0 is implied by the
            # overall "box" constraints on all the variables)
            for C_size in range(1, comb(v, k)+1):
                # Maximum number of cliques we might hit,
                # _assuming that only v vertices are present_
                # (If we pick a vertex randomly, and miss all of the
                # cliques, we get to "re-roll".)
                num_cliques_hitting_vertex = comb(v-1, k-1)
                # Bounds on number of cliques zeroed.
                # The min is at least 1 (assuming we can "re-roll"), and no
                # more than the difference between the number of cliques in the
                # larger set, and the number left over.
                min_zeroed = max(1, C_size - comb(v-1, k))
                # The max is limited by the current number of vertices (and how
                # many cliques hit one vertex), and the number of cliques.
                max_zeroed = min(num_cliques_hitting_vertex, C_size)
                # the range of possible number of cliques zeroed ...
                B_size = np.arange(min_zeroed, max_zeroed+1)
                # ... and the number left over
                A_size = C_size - B_size
                # the probability of some number of cliques being hit
                # (again, assuming that only v vertices are "in use" by
                # the hyperedges)
                def p_zonk(x):
                    return hyperg_frac(comb(v, k),
                        C_size,
                        num_cliques_hitting_vertex,
                        x)
                # the probability of at least one clique being hit
                # (if this happens, we get to "re-roll", so we assume it doesn't)
                # Note that since this we assume at least one clique is hit,
                # we normalize by this probability.
                p_at_least_one_hit = 1 - p_zonk(0)
                # ??? does this give the same answer?
                # p_at_least_one_hit = np.array([p_zonk(z1) for z1 in z]).sum()

                # The constraint on the rank of the v-vertex set in C (which is "hit")
                # is the expected rank of what's left in A, after zonking,
                # plus one (since we must have zonked at least one NAND gate)
                A = [((v, C_size), 1.)]
                A += [((v-1, C_size-j), -p_zonk(j) / p_at_least_one_hit)
                    for j in B_size]
                # pdb.set_trace()
                # We presumably "hit" a clique; depending on the basis, this
                # may mean that we also "knocked out" a gate.
                self.lp.add_constraint(A, '>', self.basis.zonked_gates())

    def add_upper_bound(self):
        """Adds upper bound.

        Here, if we have circuits for sets A and B, and let
        C = A union B, then we can construct a circuit for set C
        by ORing those two circuits together.
        """
        # loop through the number of vertices
        for v in range(self.k, self.n+1):
            # loop through number of cliques in that graph
            for C_size in range(2, comb(v, self.k)+1):
                for A_size in range(1, C_size):
                    B_size = C_size - A_size
                    # ??? on theory that it's a problem if these are the same
                    if A_size == B_size:
                        continue
                    assert(A_size >= 1)
                    assert(B_size >= 1)
                    # note that this is an _upper_ bound;
                    # also, note that three two-input NAND gates
                    # are needed to implement OR
                    self.lp.add_constraint([
                        # This is intended to make some of these
                        # constraints less linearly dependent.
                        # (Since all the variables are >= 0, this
                        # should just result in a slightly slacker
                        # bound.)
                        ((v, C_size), 1.), # np.random.uniform(1. - 0.1, 1.)),
                        ((v, A_size), -1.),
                        ((v, B_size), -1.)],
                        # XXX adding some slack here, to
                        # hopefully prevent ill-conditioned-ness
                        '<', 3.)

    def get_bounds_deprecated(self):
        """Gets bounds, for each number of vertices.

        This is the bound at each 'level'.

        ??? is there a more efficient way to compute this?
        """
        def get_bound_at_level(i):
            # this is the number we're minimizing
            x = (self.n, i)
            r = self.lp.solve(x, bounds=(0, np.inf))
            return r[x]
        # get bounds for all the vertices, and any number of cliques
        b = np.array([get_bound_at_level(c)
            for c in range(comb(self.n, self.k)+1)])
        return b

    def get_all_bounds(self):
        """Gets bounds, minimizing each variable individually."""
        def get_bound_for_variable(x):
            (num_vertices, num_cliques) = x
            # print(f"solving: n={num_vertices} num_cliques={num_cliques}")
            r = self.lp.solve(x)
            return (num_vertices, num_cliques, np.round(r, 2))
        # we only get bounds for "expected number of gates"
        bounds = [get_bound_for_variable(x)
            for x in self.expected_num_gates_vars]
        return pandas.DataFrame(bounds,
            columns = ['Num. vertices', 'Num. cliques', 'Min. gates'])

def get_bounds(n, k, max_gates, constraints_label,
        use_counting_bound, use_vertex_zeroing, use_upper_bound):
    """Gets bounds with some set of constraints.

    n, k, max_gates: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_vertex_zeroing, use_upper_bound: whether to use
        each of these constraints
    """  
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, max_gates={max_gates}, label={constraints_label}]\n')
    bound = LpVertexZeroing(n, k, max_gates)
    bound.add_level_constraints()
    if use_counting_bound:
        bound.add_counting_bound()
    if use_vertex_zeroing:
        bound.add_vertex_zeroing_constraints()
    if use_upper_bound:
        bound.add_upper_bound()
    b = bound.get_all_bounds()
    b['Constraints'] = constraints_label
    return b.iloc[:,[3,0,1,2]]

if __name__ == '__main__':
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    max_gates = int(sys.argv[3])
    # output_file = sys.argv[3]

    bounds = pandas.concat([
        get_bounds(n, k, max_gates, 'Counting', True, False, False),
        get_bounds(n, k, max_gates, 'Zeroing', False, True, False),
        get_bounds(n, k, max_gates, 'Counting and zeroing', True, True, False),
#        get_bounds(n, k, max_gates, 'All', True, True, True)
    ])
    bounds.to_csv(sys.stdout, index=False)

