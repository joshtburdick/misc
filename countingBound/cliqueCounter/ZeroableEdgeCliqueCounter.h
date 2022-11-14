#pragma once

#include <cstdlib>
#include <map>
#include <random>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "utils.h"
#include "SubsetIterator.h"
#include "CliquesByEdge.h"

using namespace std;
using namespace boost;

/** Type of bit vector implementation. */
typedef dynamic_bitset<> bits;

/** Counts (hyper)edges in a (hyper)graph, which would or wouldn't
    be "hit" by zeroing out 2-edges ("normal edges"). */
class ZeroableEdgeCliqueCounter {

public:
    /** Constructor.
        n: number of vertices in the graph
        r: size of each (hyper)edge */
    ZeroableEdgeCliqueCounter(int n, int r, CliquesByEdge & cbe);

    /** Samples a random hypergraph, then counts how many 2-edges
        could be zeroed out without hitting it (and how many
        hyperedges would be hit by those 2-edges).

        After running this, e_, zeroable_, and num_zeroable_edges_
        will all be valid. */
    void sample();

    /** Samples some graphs, and prints counts.
        numSamples: how many samples
        Side effects: prints counts to standard output */
    void printCounts(int numSamples);

private:
    /** Which cliques cover each edge. */
    CliquesByEdge & cbe_;

    /** Hyperedges in the currently-sampled hypergraph. */
    bits e_;

    /** Hyperedges which could be zeroed out (without hitting e_). */
    bits zeroable_;

    /** Number of edges which could be zeroed out. */
    int num_zeroable_edges_;
};

