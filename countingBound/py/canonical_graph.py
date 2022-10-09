# Utilities for finding canonical (hyper-)graphs, and
# related things.

import itertools

import more_itertools
import scipy.sparse.csgraph

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
        self.canonical_map = self.get_canonical_map()
        # get set of all canonical graphs
        self.canonical_graphs = frozenset(self.canonical_map.values())
        # number the canonical graphs
        self.numbering = dict(zip(self.canonical_graphs, itertools.count()))
        # get zeroing relation
        self.Z = self.get_zeroing_relation()
        # get number of graphs equivalent to each hypergraph
        self.set_size = {g: 0 for g in self.canonical_graphs}
        for k in self.canonical_map.keys():
            self.set_size[ self.canonical_map[k] ] += 1
 
    def get_canonical_map(self):
        """Constructs a dict of 'canonical hypergraphs'.

        Vertices are numbered from 0 to n-1, inclusive.
        Hypergraphs are represented as frozensets of frozensets
            of vertices.
        This is not intended to be efficient enough to work for very
            large sets of hypergraphs.
        Returns: a dict, mapping each possible hypergraph g to
            a canonical hypergraph isomorphic to g.
        Side effects: stores list of all cliques in self.all_cliques
        """
        n = self.n
        k = self.k
        # get the list of all cliques
        self.all_cliques = [frozenset(c) for c in itertools.combinations(range(n), k)]
        # This will be the mapping from a graph to a canonical (isomorphic) graph.
        # We maintain the invariant that whenever we add a graph G to this, we
        # add everything isomorphic to G as well.
        m = dict()
        # loop through the hypergraphs
        for g0 in more_itertools.powerset(self.all_cliques):
            g = frozenset(g0)
            # has g been added yet?
            if g not in m:        
                # loop through hypergraphs isomorphic to g
                for g1 in self.isomorphic_hypergraphs(g):
                    # add g1, if it hasn't been added yet
                    if g1 not in m: m[g1] = g
        # super-basic check of correctness
        assert(len(m) == 2**len(self.all_cliques))
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
        for edge in itertools.combinations(range(self.n), 2):
            for A in self.canonical_graphs:
                B = cliques_left_after_zeroing(A, frozenset(edge))
                # require that B < A (that is, the zeroing hit a clique)
                if B < A:
                    # also, use the "canonical" version of B
                    Z.add((A, self.canonical_map[B]))
        return Z

    def get_num_sets_above(self):
        """For each set, gets the number of sets above it in Z."""
        # construct sparse matrix representing Z
        n = len(self.canonical_graphs)
        M = np.zeros([n,n])
        for (A, B) in self.Z:
            M[ self.numbering[B], self.numbering[A] ] = 1
        M = csr_matrix(M)
        # get all ancestors by computing transitive closure
        dist = floyd_warshall(csgraph=M)
        pdb.set_trace()
        # for j in range(n):
        #     s = 0

