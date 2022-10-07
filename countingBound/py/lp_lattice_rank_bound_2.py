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
        self.graph_info = canonical_graph.CanonicalGraphs(n, k)
        pdb.set_trace()
        self.lp = LP_Helper(self.graph_info.canonical_graphs)







if __name__ == '__main__':
    lrb = LatticeRankBound2(4,3)   # start small




