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
        r: size of each (hyper)edge
        cbe: list of hypercliques hit by each hyperedge.
            This is separate */
    ZeroableEdgeCliqueCounter(int n, int r, CliquesByEdge & cbe);

    /** Samples a random hypergraph.

        Note that this will sample an exact number of hyperedges. Another
            option would be to simply pick hyperedges with some probability.
        num_hyperedges: the number of hyperedges the graph should have
        Side effects: sets e_. */
    void sample(int num_hyperedges);

    /** Counts how many 2-edges could be zeroed out without hitting e_
        (and how many hyperedges would be hit by those 2-edges).

        After running this, zeroable_ and num_zeroable_edges_
        will all be valid. */
    void count_zeroable_edges();

    /** Samples some graphs, and prints counts.
        num_hyperedges: how many hyperedges should be present
        num_samples: how many samples
        Side effects: prints counts to standard output */
    void printCounts(int num_hyperedges, int num_samples);

private:
    /** Number of vertices. */
    int n_;

    /** Size of each (hyper)edge. */
    int r_;

    /** Which cliques cover each edge. */
    CliquesByEdge & cbe_;

    /** Hyperedges in the currently-sampled hypergraph. */
    bits e_;

    /** Hyperedges which could be zeroed out (without hitting e_). */
    bits zeroable_;

    /** Number of edges which could be zeroed out. */
    int num_zeroable_edges_;

    /** Random number generator. */
    std::mt19937 gen_;
};

