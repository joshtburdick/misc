
#include <cstdlib>
#include <iostream>

#include "ZeroableEdgeCliqueCounter.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 4) {
        cerr << "usage: zeroable_edges n r numSamples" << endl;
        cerr << "n: number of vertices" << endl;
        cerr << "r: number of vertices per (hyper)edge" << endl;
        cerr << "numSamples: number of samples" << endl;
        return 1;
    }

    int n = atoi(argv[1]);
    int r = atoi(argv[2]);
    int num_samples = atoi(argv[3]);

    if (!(n > r)) {
        cerr << "error: need n > r, but n = " << n << " and r = " << r << endl;
        return 1;
    }
    CliquesByEdge cbe(n, r);
    cout << "constructed CliquesByEdge" << endl;
    ZeroableEdgeCliqueCounter zeroable_edges(n, r, cbe);
    cout << "constructed ZeroableEdgeCliqueCounter" << endl;
    zeroable_edges.printCounts(num_samples);
    return 0;
}

