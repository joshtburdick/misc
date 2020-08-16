/** Counts cliques, and how many edges they cover,
    in random hypergraphs. */

#include <cstdlib>
#include <iostream>

#include "CliqueCoverageCounter.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 5) {
        cerr << "usage: coverage n r maxCliqueSize numSamples" << endl;
        cerr << "n: number of vertices" << endl;
        cerr << "r: number of vertices per (hyper)edge" << endl;
        cerr << "maxCliqueSize: only look for cliques up to this size" << endl;
        cerr << "  (can be much smaller than n)" << endl;
        cerr << "numSamples: number of (pairs of) samples" << endl;
        cerr << "  (as this samples graphs and their complements" << endl;
        return 1;
    }

    int n = atoi(argv[1]);
    int r = atoi(argv[2]);
    int maxCliqueSize = atoi(argv[3]);
    int numSamples = atoi(argv[4]);

    if (!(n > r)) {
        cerr << "error: need n > r" << endl;
        return 1;
    }
    if (!(r <= maxCliqueSize && maxCliqueSize <= n)) {
        cerr << "error: need r <= maxCliqueSize <= n" << endl;
        return 1;
    }

    CliqueCoverageCounter counter(n, r, maxCliqueSize);
    for(int s = 0; s <= numSamples; ++s)
        counter.printCoverageCounts();
    return 0;
}

