#include "MaxCliqueCounter.h"

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

/** Constructor. */
MaxCliqueCounter::MaxCliqueCounter(int n, int r) :
        cliques(n, r, n), count(n+1) {
}

/** Samples some random graphs, and prints counts. */
void MaxCliqueCounter::printMaxCounts(int numSamples) {
    // clear the counts
    for(int s = 0; s < cliques.getNumVertices(); s++)
        count[s] = 0;
    // sample a uniformly random graph, find size of the
    // largest (hyper)clique, and add to the counts
    for(int iter = 0; iter < numSamples; ++iter) {
        bits g = randomBitset(cliques.getNumEdges());
        int maxCliqueSize = cliques.getMaxCliqueSize(g);
        ++count[ maxCliqueSize ];
    }

    // print the counts
    for(int s = 0; s < cliques.getNumVertices(); s++)
        cout << count[s] << " ";
    cout << endl;
}

/** Creates a random bitset of some size. */
bits MaxCliqueCounter::randomBitset(int numBits) {
    bits b;
    // first, get enough 64-bit blocks
    for(int i = 0; i < (numBits/64)+1; i++) {
        // get a 64-bit random number
        uint64_t r = random() ^ (random() << 32);
        b.append(r);
    }
    // then truncate to the required size
    b.resize(numBits);
    return b;
}

