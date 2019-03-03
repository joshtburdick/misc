
#include <cstdlib>
#include <iostream>

#include "MaxCliqueCounter.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 4) {
        cerr << "usage: count n r numSamples" << endl;
        cerr << "n: number of vertices" << endl;
        cerr << "r: number of vertices per (hyper)edge" << endl;
        cerr << "numSamples: number of samples, except:" << endl;
        cerr << "  -1: loop through all possible hypergraphs"
            << " (not yet implemented)" << endl;
        cerr << "  -2: dump hypercliques (for testing), and exit" << endl;
        return 1;
    }

    int n = atoi(argv[1]);
    int r = atoi(argv[2]);
    int numSamples = atoi(argv[3]);

    if (!(n > r)) {
        cerr << "error: need n > r, but n = " << n << " and r = " << r << endl;
        return 1;
    }
    if (numSamples == -1) {
        cerr << "error: looping through all hypergraphs "
            << "not yet implemented" << endl;
    }
    // dump the cliques
    if (numSamples == -2) {
        CliqueList cliques(n, r, n);
        cliques.printCliques();
    }
    if (numSamples > 0) {
        MaxCliqueCounter counter(n, r);
        while (1)
            counter.printMaxCounts(numSamples);
    }
    return 0;
}

