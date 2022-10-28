#pragma once

#include <cstdlib>
#include <map>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "SubsetIterator.h"

using namespace std;
using namespace boost;

/** Type of bit vector implementation. */
typedef dynamic_bitset<> bits;
using namespace std;

/** Counts (hyper)cliques, after input edges are zeroed out. */
class ZeroedEdgeCliqueCounter {
public:
    /** Constructor.
        n: number of vertices in the graph
        r: size of each (hyper)edge */
    MaxCliqueCounter(int n, int r);

    /** Samples some graphs, and prints counts.
        numSamples: how many samples
        Side effects: prints counts to standard output */
    void printMaxCounts(int numSamples);

private:
    /** List of cliques to check for. */
    CliqueList cliques;

    /** The counts of cliques of each size. */
    vector<int> count;        

    /** Creates a random bitset of some size. */
    bits randomBitset(int numBits);
};

