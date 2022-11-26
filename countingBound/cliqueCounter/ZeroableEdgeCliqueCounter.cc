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
    zeroable_(cbe_.hyperedges_.size()), gen_() {

    // initialize random number source
    // ??? is this the right way to do this?
    std::random_device rd;
    gen_.seed(rd());
}

void ZeroableEdgeCliqueCounter::sample(int num_hyperedges) {
    // first, initialize bits, based on the number of hyperedges
    // ??? is this the simplest way to make p floating-point?
    double p = (1.0 * num_hyperedges) / (1.0 * e_.size());
    std::bernoulli_distribution random_bit(p);
    for(int i=0; i<e_.size(); ++i) {
        bool bit = random_bit(gen_);
        e_[i] = (bit ? 1 : 0);
    }
    // then, tweak this to match num_hyperedges
    std::uniform_int_distribution<> random_edge(0, e_.size()-1);
    int num_bits_set = e_.count();
    while (num_bits_set != num_hyperedges) {
        if (num_bits_set < num_hyperedges)
            e_.set(random_edge(gen_));
        if (num_bits_set > num_hyperedges)
            e_.reset(random_edge(gen_));
        num_bits_set = e_.count();
    }
    // at this point, e_ should have num_hyperedges bits set
}

void ZeroableEdgeCliqueCounter::count_zeroable_edges() {
    // clear mask of hyperedges which would be "hit" by those
    zeroable_.reset();
    // clear count of which 2-edges could be zeroed
    num_zeroable_edges_ = 0;
    // loop through the 2-edges
    // cout << "e_ = " << e_ << endl;
    for(int i=0; i<cbe_.edge_mask_.size(); ++i) {
        bits & z = cbe_.edge_mask_[i];
        // cout << "z = " << z << endl;
        // if none of these edges are in the graph
        if ((e_ & z).none()) {
            // record that this 2-edge could be zeroed
            ++num_zeroable_edges_;
            zeroable_ |= z;
        }
    }
}

void ZeroableEdgeCliqueCounter::printCounts(int num_hyperedges, int num_samples) {
    for(int iter=0; iter<num_samples; ++iter) {
        sample(num_hyperedges);
        count_zeroable_edges();
        cout << e_.count() << ','
            << zeroable_.count() << ','
            << num_zeroable_edges_ << endl;
    }
}

