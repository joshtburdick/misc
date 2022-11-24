#include "ZeroableEdgeCliqueCounter.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "utils.h"
#include "SubsetIterator.h"

using namespace std;
using namespace boost;

/** Type of bit vector implementation. */
typedef dynamic_bitset<> bits;

ZeroableEdgeCliqueCounter::ZeroableEdgeCliqueCounter(
        int n, int r, CliquesByEdge & cbe) :
    n_(n), r_(r), cbe_(cbe), e_(1), zeroable_(1) {

    cout << "created ZECC" << endl;
}

void ZeroableEdgeCliqueCounter::sample() {
    // sample a random hypergraph
    e_ = randomBitset(e_.size());
    cout << "sampled hypergraph" << endl;
    // clear mask of hyperedges which would be "hit" by those
    zeroable_.reset();
    // clear count of which 2-edges could be zeroed
    num_zeroable_edges_ = 0;
    // loop through the 2-edges
    for(int i=0; i<cbe_.edge_mask_.size(); ++i) {
        bits & z = cbe_.edge_mask_[i];
        // if these edges aren't in the graph
        if ((e_ & z).none()) {
            // record that this 2-edge could be zeroed
            ++num_zeroable_edges_;
            zeroable_ |= z;
        }
    }
}

void ZeroableEdgeCliqueCounter::printCounts(int numSamples) {
    for(int iter=0; iter<numSamples; ++iter) {
        cout << e_.count() << ','
            << zeroable_.count() << ','
            << num_zeroable_edges_ << endl;
    }
}

