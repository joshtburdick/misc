#pragma once

#include <vector>

#include <cstdlib>

#include "CliqueList.h"

using namespace std;

/** Counts how often the largest (hyper)clique was a given size. */
class CliqueCoverageCounter {
public:
    /** Constructor.
        n: number of vertices in the graph
        r: size of each (hyper)edge
	maxCliqueSize: the largest size of clique to check for */
    CliqueCoverageCounter(int n, int r, int maxCliqueSize);

    /** Samples a graph, and prints counts of cliques of various
        sizes (and how many edges they cover). Then, prints counts
        for the graph's complement.
        Side effects: prints counts to standard output */
    void printCoverageCounts();

private:
    /** List of cliques to check for. */
    CliqueList cliques;

    /** The counts of cliques of each size. */
    vector<int> cliqueCount;        

    /** The counts of edges covered by each size of clique. */
    vector<int> coveredEdgeCount;        

    /** Creates a random bitset of some size. */
    bits randomBitset(int numBits);

    /** Prints coverage counts for one graph. */
    void printCoverageCounts1(bits & g);
};

