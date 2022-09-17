# Utilities for finding canonical (hyper-)graphs.

import itertools

import more_itertools


def isomorphic_hypergraphs(n, g):
    """Enumerates all hypergraphs isomorphic to a hypergraph."""
    def renumber_edge(pi, edge):
        return frozenset([pi[v] for v in edge])
    def renumber_graph(pi, g):
        return frozenset([renumber_edge(pi, edge) for edge in g])
    return [renumber_graph(pi, g)
        for pi in itertools.permutations(range(n))]

def canonical_hypergraph_map(n, k):
    """Constructs a dict of 'canonical hypergraphs'.

    Vertices are numbered from 0 to n-1, inclusive.
    Hypergraphs are represented as frozensets of frozensets
        of vertices.
    This is not intended to be efficient enough to work for very
        large sets of hypergraphs.
    n: number of vertices in the graph
    k: number of vertices in a clique
    Returns: a dict, mapping each possible hypergraph g to
        a canonical hypergraph isomorphic to g.
    """
    # get the list of all cliques
    all_cliques = [frozenset(c) for c in itertools.combinations(range(n), k)]
    # This will be the mapping from a graph to a canonical (isomorphic) graph.
    # We maintain the invariant that whenever we add a graph G to this, we
    # add everything isomorphic to G as well.
    m = dict()
    # loop through the hypergraphs
    for g in more_itertools.powerset(all_cliques):
        # has g been added yet?
        if g not in m:        
            # loop through hypergraphs isomorphic to g
            for g1 in isomorphic_hypergraphs(n, g):
                # add g1, if it hasn't been added yet
                if g1 not in m: m[g1] = g
    # super-basic check of correctness
    assert(len(m) == 2**len(all_cliques))
    return m




