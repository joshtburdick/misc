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
    n_(n), r_(r), cbe_(cbe),
    e_(cbe_.hyperedges_.size()),
    zeroable_(cbe_.hyperedges_.size()) {

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
    cout << "e_ = " << e_ << endl;
    for(int i=0; i<cbe_.edge_mask_.size(); ++i) {
        bits & z = cbe_.edge_mask_[i];
        cout << "z = " << z << endl;
        // if none of these edges are in the graph
        if ((e_ & z).none()) {
            // record that this 2-edge could be zeroed
            ++num_zeroable_edges_;
            zeroable_ |= z;
        }
    }
}

void ZeroableEdgeCliqueCounter::printCounts(int numSamples) {
    for(int iter=0; iter<numSamples; ++iter) {
        sample();
        cout << e_.count() << ','
            << zeroable_.count() << ','
            << num_zeroable_edges_ << endl;
    }
}

