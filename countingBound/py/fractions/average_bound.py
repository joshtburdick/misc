#!/usr/bin/env python3
# Attempt at average bounds on circuit size, for large random sets of cliques.
# FIXME
# - in level constraints, the expected value variables should include all
#   functions with <= that number of vertices.
# - use argparse
# - compute "canonical" solution (first finding the minimum for the objective
#   function, then finding a solution with the minimum total expected number of
#   gates, which has the same value for the objective function).

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
import pulp_helper

sys.path.append("..")

import hypergraph_counter

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

        This will have several groups of variables:
        - Tuples of the form ("g", v, c, gates) where:
            - `v` is as above
            - `c` is as above
            - `gates` is the number of gates
            - Each variable will be the number of sets of cliques
              (with `c` cliques and _exactly_ `v` vertices), with smallest
              circuit having exactly `gates` gates.
        - Tuples of the form ("v", v, c) where:
            - `v` is the number of vertices, with k <= v <= n
            - `c` is the number of cliques, with 0 <= c <= {n choose v}
            - Each variable will be the expected number of gates in
            the sets of cliques of size `c`, using _exactly_ `v` vertices.
        - Tuples of the form ("u", v, c) which are like "v", but for
          the case where the sets of cliques use _up to_ `v` vertices.

        n: number of vertices in the graph
        k: number of vertices in a clique (>= 3)
        max_gates: maximum number of gates to include
        """
        self.n = n
        self.k = k
        if k < 3:
            raise ValueError('k must be >= 3')
        self.max_gates = max_gates
        # number of possible cliques (in the k-uniform hypergraph on n vertices)
        self.num_possible_cliques = comb(n, k)
        
        # counts of hypergraphs with exactly some number of vertices
        self.hc = hypergraph_counter.HypergraphCounter(self.n, self.k)
        self.hypergraph_counts = self.hc.count_hypergraphs_exact_vertices()

        # Variables for expected number of gates; initially just the 0-clique case.
        self.expected_num_gates_vars = []
        for v in range(k-1, n+1):
            self.expected_num_gates_vars += [("v", v, 0), ("u", v, 0)]

        # [('v', k-1, 0), ('u', k-1, 0)]
        # Variables for counts of numbers of functions (for which we don't
        # include the 0-clique case).
        self.num_gates_dist_vars = []
        for v in range(k, n+1):
            for c in range(comb(v, k) + 1):
                self.expected_num_gates_vars += [("v", v, c), ("u", v, c)]
                self.num_gates_dist_vars += [("g", v, c, g) for g in range(1, max_gates+1)]
        # presumably, since these counts start at 1 gste, we don't need to constrain
        # the expected number of gates to be >= 1
        
        # wrapper for LP solver
        self.lp = pulp_helper.PuLP_Helper(
            self.expected_num_gates_vars + self.num_gates_dist_vars)

        # self.basis = gate_basis.UnboundedFanInNandBasis()
        self.basis = gate_basis.TwoInputNandBasis()

    def add_averaging_constraints(self):
        """Adds averaging constraints.

        These include:
        - constraining the total number of functions in each group
          (from `self.hypergraph_counts`)
        - connecting `v`, the expected number of gates in each group 
          and `g`, the number of functions with a given number of gates
        - connecting `v` and `u` (since `u` averages together several
          levels of `v`)
        """
        # loop through the groups of functions
        for v in range(self.k, self.n+1):
            num_possible_cliques = comb(v, self.k)
            for c in range(num_possible_cliques+1):
                num_functions = self.hypergraph_counts[v][c]
                # add constraint that these sum to the number of functions
                self.lp.add_constraint(
                    [(("g", v, c, i), 1) for i in range(1, self.max_gates+1)],
                    '=', num_functions)
                # add constraint defining expected number of gates
                A = [(("g", v, c, i), i)
                    for i in range(1, self.max_gates+1)]
                self.lp.add_constraint(A + [(("u", v, c), -num_functions)],
                    '=', 0)

        # add constraints connecting "u" and "v"
        for max_v in range(self.k, self.n+1):
            num_possible_cliques = comb(max_v, self.k)
            for c in range(num_possible_cliques + 1):
                A = []
                total_functions = 0
                for v in range(self.k, max_v+1):
                    if c < self.hypergraph_counts[v].shape[0]:
                        num_functions = self.hypergraph_counts[v][c]
                        A.append((("v", v, c), num_functions))
                        total_functions += num_functions
                self.lp.add_constraint([(("u", max_v, c), -total_functions)] + A, '=', 0)

        # add constraints for sets of zero cliques
        for v in range(self.k-1, self.n+1):
            self.lp.add_constraint([(("u", v, 0), 1)], "=", 1)
            self.lp.add_constraint([(("v", v, 0), 1)], "=", 1)

    def add_counting_bound(self):
        """Adds counting bounds, for a given number of gates.

        This is an upper bound on the number of functions with
        some number of gates. (Hopefully, this will force some functions
        to need a larger number of gates, giving a lower bound.)

        Note that many sets of cliques will be the same, up to isomorphism;
        we count these as separate. (This is because, for instance, if n=6,
        if we zero out one vertex, the remaining cliques count as distinct circuits,
        even though they're the same, up to edge labelling.)
        """
        num_functions = self.basis.num_functions(self.num_possible_cliques, self.max_gates)
        # note that this doesn't include the empty set of cliques
        for g in range(1, self.max_gates+1):
            self.lp.add_constraint([(("g", v, c, g), 1)
                for (t, v, c) in self.expected_num_gates_vars
                if t == "u" and v >= self.k],
                '<=', num_functions[g]-1)

    def add_vertex_zeroing_constraints(self, use_lower_bound, use_upper_bound):
        """Adds constraints from zeroing out vertices.

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
            for C_size in range(1, comb(v, self.k)+1):
                # Maximum number of cliques we might hit,
                # _assuming that only v vertices are present_.
                # (If we pick a vertex randomly, and miss all of the
                # cliques, we get to "re-roll".)
                num_cliques_hitting_vertex = comb(v-1, self.k-1)
                # Bounds on number of cliques zeroed.
                # The min is at least 1 (assuming we can "re-roll"), and no
                # more than the difference between the number of cliques in the
                # larger set, and the number left over.
                min_zeroed = max(1, C_size - comb(v-1, self.k))
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
                p_hit = np.array([hyperg_frac(
                    comb(v, self.k),
                    C_size,
                    num_cliques_hitting_vertex,
                    h)
                    for h in B_size])
                # The above calculation omits the possibility of missing all of the cliques.
                # We argue that we can ignore this possibility (since for a non-empty set of
                # cliques, there's always _some_ vertex which will hit >= 1 clique).)
                # Therefore, we normalize this (by dividing by the sum).
                p_hit /= p_hit.sum()
                # These coefficients are the difference between the expected number of gates
                # in A (before zeroing out a vertex) and C (after zeroing out a vertex).
                A = [(("u", v, C_size), 1)]
                A += [(("u", v-1, A_size[i]), -p_hit[i]) for i in range(B_size.size)]
                # pdb.set_trace()
                # Lower bound: we "hit" a clique, and so we must have zonked at least
                # one NAND gate. (For other bases, this might not be guaranteed...)
                if use_lower_bound:
                    self.lp.add_constraint(A, ">=", self.basis.zonked_gates())
                # Upper bound: The number of gates "zonked" (in B) is no more than
                # the number of cliques "hit" (in B), since we could have implemented
                # C by taking the circuit for A, and "patching" it to include the cliques in B,
                # (Note that this only holds for the unbounded-fan-in NAND gate basis.)
                upper_bound = (p_hit * B_size).sum()
                print(f"v={v}, C_size={C_size}, B_size={B_size}, p_hit={p_hit}, upper_bound={upper_bound}")
                if use_upper_bound:
                    self.lp.add_constraint(A, "<=", upper_bound)

    def get_all_bounds(self):
        """Gets bounds for each possible number of cliques.

        This is the bounds for each possible number of cliques,
        in the scenario that the number of gates for CLIQUE (or functions
        detecting many cliques) is minimized.
        """
        # solve, minimizing number of gates for CLIQUE
        r = self.lp.solve(("u", self.n, self.num_possible_cliques))
        if not r:
            return None
        # for now, we only get bounds for "expected number of gates"
        # for each number of cliques
        n_cliques = range(self.num_possible_cliques+1)
        bounds = [r[("u", self.n, c)] for c in n_cliques]
        return pandas.DataFrame({
                'Num. vertices': self.n,
                'Num. cliques': n_cliques,
                'Min. gates': bounds})
        # FIXME return value of objective function?

def get_bounds(n, k, max_gates, constraints_label,
        use_counting_bound, use_zeroing_lower_bound, use_zeroing_upper_bound):
    """Gets bounds with some set of constraints.

    n, k, max_gates: problem size
    constraints_label: label to use for this set of constraints
    use_counting_bound, use_zeroing_lower_bound, use_zeroing_upper_bound:
        whether to use each of these constraints
    Returns:
        A pandas DataFrame with columns ['Num. vertices', 'Num. cliques', 'Min. gates']
    """  
    # ??? track resource usage?
    sys.stderr.write(f'[bounding with n={n}, k={k}, max_gates={max_gates}, label={constraints_label}]\n')
    bound = LpVertexZeroing(n, k, max_gates)
    bound.add_averaging_constraints()
    if use_counting_bound:
        bound.add_counting_bound()
    bound.add_vertex_zeroing_constraints(use_zeroing_lower_bound, use_zeroing_upper_bound)

    b = bound.get_all_bounds()
    b['Constraints'] = constraints_label
    return b.iloc[:,[3,0,1,2]]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int, help='number of vertices')
    parser.add_argument('k', type=int, help='size of cliques')
    parser.add_argument('max_gates', type=int, help='maximum number of gates')
    parser.add_argument('band_width', type=int, help='width of highest possible number of cliques to minimize')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args() 
    # output_file = sys.argv[3]

    bounds = pandas.concat([
        get_bounds(args.n, args.k, args.max_gates, 'Counting', True, False, False),
        get_bounds(args.n, args.k, args.max_gates, 'Counting, zeroing lower', True, True, False),
        get_bounds(args.n, args.k, args.max_gates, 'Counting, zeroing upper', True, False, True),
        get_bounds(args.n, args.k, args.max_gates, 'Counting, zeroing lower and upper', True, True, True),
    ])
    bounds.to_csv(sys.stdout, index=False)
