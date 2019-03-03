#pragma once

#include <map>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "SubsetIterator.h"

using namespace std;
using namespace boost;

/** Type of bit vector implementation. */
typedef dynamic_bitset<> bits;

/** List of cliques occurring in a larger graph, able to count
    the size of the largest clique in an arbitrary graph. */
class CliqueList {
public:
    /** The number of vertices. */
    int n_;

    /** The size of the edges (as this works for hypergraphs). */
    int r_;

    /** The largest size of cliques to include. */
    int max_clique_size_;

    /** Constructor. Note that this isn't particularly fast (but
        doesn't have to be). */
    CliqueList(int n, int r, int maxCliqueSize);

    /** Adds counts of which cliques are present in a graph.
        Should be multithreaded.
        Possibly not used, as the bookkeeping seems easier if I
        just track the largest clique. (Clearly, if there's a
        9-vertex clique, there are many 8-cliques, and even
        more 7-vertex cliques, etc.) */
    void addToCounts(vector<int> count, const bits & g);

    /** Gets size of the largest clique in a graph.
        Should be multithreaded.
        g: vertices of a hypergraph
        Returns: number of vertices in the largest clique */
    int getMaxCliqueSize(const bits & g);

    /** Gets vertices covered by hypergraphs.
        g: vertices of a hypergraph
        k: size of the hyperclique to check for
        coveredEdges: ref. to bits, which will be set with
            the edges which are covered
        Side effects: fills in coveredEdges
        Returns: number of k-vertex cliques found */
    int getCoveredEdges(const bits & g, int k, bits & coveredEdges);

    /** Dumps the list of cliques to standard output (for debugging). */
    void printCliques();

    /** Gets the number of vertices. */
    int getNumVertices() {
        return n_;
    }

    /** Gets the number of potential edges. */
    int getNumEdges() {
        return edgeIndex_.size();
    }

private:

    /** Mapping from edge to index in the bit vectors. */
    map<vector<int>, int> edgeIndex_;

    /** The cliques which we're testing for.
        This is indexed first by the number of vertices in the clique;
        each element of this is a vector of the cliques. */
    vector<vector<bits> > clique_;

    /** Initialize the edge numbering. */
    void initEdgeNumbering();

    /** Initialize the cliques with k vertices. */
    void addCliques(int k);

    /** Gets the bitmask corresponding to a particular clique. */
    bits getBitmask(vector<int> & edge);

    /** Prints the edges in one graph (by translating the bits which
        are set to vertex numbers). */
    void printGraph(bits &e);
};

