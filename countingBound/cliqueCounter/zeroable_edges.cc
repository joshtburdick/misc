
#include <cstdlib>
#include <iostream>

#include "ZeroableEdgeCliqueCounter.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 5) {
        cerr << "usage: zeroable_edges n r edges numSamples" << endl;
        cerr << "n: number of vertices" << endl;
        cerr << "r: number of vertices per (hyper)edge" << endl;
        cerr << "numEdges: number of (hyper)edges" << endl;
        cerr << "numSamples: number of samples" << endl;
        return 1;
    }

    int n = atoi(argv[1]);
    int r = atoi(argv[2]);
    int num_hyperedges = atoi(argv[3]);
    int num_samples = atoi(argv[4]);

    if (!(n > r)) {
        cerr << "error: need n > r, but n = " << n << " and r = " << r << endl;
        return 1;
    }
    CliquesByEdge cbe(n, r);
    if (!(num_hyperedges <= cbe.hyperedges_.size())) {
        cerr << "error: cannot sample " << num_hyperedges << " hyperedges, ";
        cerr << "since there are only " << cbe.hyperedges_.size() << " possible hyperedges" << endl;
        return 1;
    }
    // cout << "constructed CliquesByEdge" << endl;
    ZeroableEdgeCliqueCounter zeroable_edges(n, r, cbe);
    // cout << "constructed ZeroableEdgeCliqueCounter" << endl;
    zeroable_edges.printCounts(num_hyperedges, num_samples);
    return 0;
}

