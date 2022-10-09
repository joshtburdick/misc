#!/usr/bin/env python3
# Bound using the lattice of canonical hypergraphs.

import pdb

import canonical_graph
import lp_helper

class LatticeRankBound2:
    """Bound, using the lattice of canonical hypergraphs.

    """

    def __init__(self, n, k):
        """Constructor gets graph info, and sets up variable names.

        n: number of vertices in the graph
        k: number of vertices in a clique
        """
        # get info about the the graphs in the lattice
        self.graph_info = canonical_graph.CanonicalGraphs(n, k)
        # set up linear program, with one variable per (class of) graph
        self.lp = lp_helper.LP_Helper(self.graph_info.canonical_graphs)

    def add_average_rank_constraint(self):
        """Adds equality constraint on average rank of all of the functions."""
        num_sets = sum(self.graph_info.set_size.values())
        A = [(S, n / num_sets)
            for (S, n) in self.graph_info.set_size.items()]
        self.lp.add_constraint(A, '=', num_sets / 2.)

    def add_zeroing_constraints(self):
        """Add constraints that zeroing sets reduces the rank."""
        for (A,B) in self.graph_info.Z:
            self.lp.add_constraint([(A, 1.), (B, -1.)], '>', 1)

    def add_higher_sets_constraints(self):
        """Add upper bound constraints, based on "higher" sets.

        For each set A, we find the number of graphs "above and including"
        that set; call that number b.
        We then know that S is "below" all of those sets. Since we also
        know |A|, we get that E[|C(A)|] <= (N - b) - |A|/2 .
        """
        z = self.graph_info.get_num_higher_sets()

    def get_all_set_bounds(self):
        """Gets bounds for all the sets.

        This solves the linear program, including all constraints which
            have been included.
        Returns: a hash, keyed by canonical set, of the lower bound
            for the rank of each set.
        """
        cliques = self.graph_info.all_cliques
        bounds = self.lp.solve(frozenset(cliques),
            bounds=(0, 2**len(cliques)-1))
        return bounds

    def get_clique_bound(self):
        """Gets the bound for finding all the cliques."""
        bounds = self.get_all_set_bounds()
        return bounds[frozenset(self.graph_info.all_cliques)]

if __name__ == '__main__':
    lrb = LatticeRankBound2(4,3)   # start small, eh
    lrb.add_average_rank_constraint()
    lrb.add_zeroing_constraints()
    lrb.add_higher_sets_constraints()
    print(str(lrb.get_clique_bound()))

