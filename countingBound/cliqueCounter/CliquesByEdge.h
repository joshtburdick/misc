#pragma once

#include <cstdlib>
#include <map>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "utils.h"
#include "SubsetIterator.h"

using namespace std;
using namespace boost;

/** Computes which hypercliques are "missed" by each edge. */
class CliquesByEdge {
public:
    /** Constructor.
        n: number of vertices in the graph
        r: size of each (hyper)edge */
    CliquesByEdge(int n, int r);

    /** The edges, ordered by their index. */
    vector<vector<int> > edges;

    /** The hyperedges, ordered by their index. */
    vector<vector<int> > hyperedge;

    /** Masks of which hyperedges _aren't_ hit by each edge.

        This is indexed by edge index, while the bits
    are indexed by hyperedge index. */
    vector<bits> edge_mask_;

    /** Number of vertices. */
    int n_;

    /** Size of each (hyper)edge. */
    int r_;

private:
};
