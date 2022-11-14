#include "CliquesByEdge.h"

#include <cstdlib>
#include <map>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "SubsetIterator.h"

using namespace std;
using namespace boost;

/** Constructor.
    n: number of vertices in the graph
    r: size of each (hyper)edge

    If this were multithreaded, it would be faster.
    However, that may not matter much.
    */
CliquesByEdge::CliquesByEdge(int n, int r) {
    n_ = n;
    r_ = r;
    // fill in list of edges
    SubsetIterator edgeIterator(n_, 2);
    do {
        edges_.push_back(edgeIterator.getSet());
    } while (edgeIterator.next());
    // similarly, fill in list of hyperedges
    SubsetIterator hyperedgeIterator(n_, r_);
    do {
        hyperedges_.push_back(hyperedgeIterator.getSet());
    } while (hyperedgeIterator.next());
    // loop through the edges
    for(int i=0; i<edges_.size(); i++) {
        vector<int> e = edges_[i];
        // create bitset of edges
        bits m(hyperedges_.size());
        // loop through the hyperedges
        for(int j=0; j<hyperedges_.size(); j++) {
            vector<int> h = hyperedges_[j];
            // does this 2-edge "hit" this hyperedge?
            if (count(h.begin(), h.end(), e[0]) &&
                    count(h.begin(), h.end(), e[1]))
                // if so, set the corresponding bit
                m[j] = 1;
        }
        // store the mask for this edge
        edge_mask_[i] = m;
    } 
}

