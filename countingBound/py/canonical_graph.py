# Utilities for finding canonical (hyper-)graphs, and
# related things.

import itertools

import more_itertools


def cliques_left_after_zeroing(clique_set, edge):
    """Finds cliques which are left after zeroing out an edge.

    clique_set: a set of cliques
    edge: an edge (as a two-element set)
    Returns: the cliques in clique_set which aren't hit by that edge.
    """
    return frozenset([c for c in clique_set if not edge < c])

class CanonicalGraphs:
    """Gets canonical hypergraphs, and related things.
    """

    def __init__(self, n, k):
        """Constructor calculates various things.

        This includes a mapping to canonical hypergraphs,
        a count of the graphs isomorphic to each of those,
        and relations related to zeroing edges.
        """
        self.n = n
        self.k = k
        # compute mapping from graphs to a canonical graph
        self.canonical_map = self.canonical_hypergraph_map()
        # get set of all canonical graphs
        self.canonical_graphs = frozenset(self.canonical_map.values())
        # get number of graphs equivalent to each hypergraph
        self.set_size = {g: 0 for g in self.canonical_graphs()}
        for k in self.canonical_map.keys():
            self.set_size[ self.canonical_map[k] ] += 1
        # get zeroing relation
        self.Z = self.get_zeroing_relation()

    def canonical_hypergraph_map(self):
        """Constructs a dict of 'canonical hypergraphs'.

        Vertices are numbered from 0 to n-1, inclusive.
        Hypergraphs are represented as frozensets of frozensets
            of vertices.
        This is not intended to be efficient enough to work for very
            large sets of hypergraphs.
        Returns: a dict, mapping each possible hypergraph g to
            a canonical hypergraph isomorphic to g.
        """
        n = self.n
        k = self.k
        # get the list of all cliques
        all_cliques = [frozenset(c) for c in itertools.combinations(range(n), k)]
        # This will be the mapping from a graph to a canonical (isomorphic) graph.
        # We maintain the invariant that whenever we add a graph G to this, we
        # add everything isomorphic to G as well.
        m = dict()
        # loop through the hypergraphs
        for g0 in more_itertools.powerset(all_cliques):
            g = frozenset(g0)
            # has g been added yet?
            if g not in m:        
                # loop through hypergraphs isomorphic to g
                for g1 in isomorphic_hypergraphs(n, g):
                    # add g1, if it hasn't been added yet
                    if g1 not in m: m[g1] = g
        # super-basic check of correctness
        assert(len(m) == 2**len(all_cliques))
        return m

    def isomorphic_hypergraphs(self, g):
        """Enumerates all hypergraphs isomorphic to a hypergraph."""
        def renumber_edge(pi, edge):
            return frozenset([pi[v] for v in edge])
        def renumber_graph(pi, g):
            return frozenset([renumber_edge(pi, edge) for edge in g])
        return [renumber_graph(pi, g)
            for pi in itertools.permutations(range(self.n))]

    def get_zeroing_relation(self):
        """Gets the relation showing the result of zeroing an edge.

        (A, B) is in the relation Z if zeroing some edge in A yields B.
        Here, A and B are canonical hypergraphs.
        """
        Z = set()
        for edge in itertools.combinations(range(self.n)):
            for A in self.canonical_graphs:
                B = cliques_left_after_zeroing(A, edge)
                Z += (A, B)





