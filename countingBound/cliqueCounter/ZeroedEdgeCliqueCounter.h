#pragma once

#include <cstdlib>
#include <map>
#include <random>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "utils.h"
#include "SubsetIterator.h"

using namespace std;
using namespace boost;

/** Type of bit vector implementation. */
typedef dynamic_bitset<> bits;

/** Counts (hyper)edges in a (hyper)graph, after 2-edges
    ("normal edges") are zeroed out. */
class ZeroedEdgeCliqueCounter {
public:
    /** Constructor.
        n: number of vertices in the graph
        r: size of each (hyper)edge */
    MaxCliqueCounter(int n, int r, CliquesByEdge cbe &) :
        n_(n), r_(r), cbe(cbe);

    /** Samples a random hypergraph (and sets the mask of
        "hyperedges not hit by an edge" to all 1's)
        Also sets the edges to a random permutation. */
    void randomize();

    /** Picks a random edge, and zeros out hyperedges
        which contain it.
        i: index of the edge to zero */
    void zeroEdge(int i);

    /** Samples a random hypergraph, then zeros out edges
        until the graph is empty. */
    void sampleGraphsDownToZero();

    /** Samples some graphs, and prints counts.
        numSamples: how many samples
        Side effects: prints counts to standard output */
    void printMaxCounts(int numSamples);

private:
    /** Number of vertices. */
    int n_;

    /** Size of each (hyper)edge. */
    int r_;

    /** Which cliques cover each edge. */
    CliquesByEdge cbe &;

    /** Hyperedges in the currently-sampled hypergraph. */
    bitset e_;

    /** Hyperedges which haven't been zeroed out yet. */
    bitset not_zeroed_;

    /** Order in which edges should be zeroed out. */
    vector<int> edge_order_;
};
