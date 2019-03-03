#include "CliqueCoverageCounter.h"

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

/** Constructor. */
CliqueCoverageCounter::CliqueCoverageCounter(int n, int r, int maxCliqueSize) :
    cliques(n, r, maxCliqueSize),
    cliqueCount(maxCliqueSize+1),
    coveredEdgeCount(maxCliqueSize+1) {
}

/** Samples a random graph, and prints counts of cliques of various sizes
    (and how many edges they cover). */
void CliqueCoverageCounter::printCoverageCounts() {
    // sample a graph, and print counts
    bits g = randomBitset(cliques.getNumEdges());
    printCoverageCounts1(g);    
    // then print counts for the complement
    g.flip();
    printCoverageCounts1(g);
}

/** Prints counts for one graph. */
void CliqueCoverageCounter::printCoverageCounts1(bits & g) {
    bits coveredEdges(g);
    for(int k = 0; k <= cliques.max_clique_size_; ++k) {
        cliqueCount[k] = cliques.getCoveredEdges(g, k, coveredEdges);
        coveredEdgeCount[k] = coveredEdges.count();
    }
    // print the counts
    for(int k = 0; k <= cliques.max_clique_size_; ++k)
        cout << cliqueCount[k] << " ";
    cout << " ";  // add a tiny bit of formatting
    for(int k = 0; k <= cliques.max_clique_size_; ++k)
        cout << coveredEdgeCount[k] << " ";
    cout << endl;
}

/** Creates a random bitset of some size. */
bits CliqueCoverageCounter::randomBitset(int numBits) {
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

