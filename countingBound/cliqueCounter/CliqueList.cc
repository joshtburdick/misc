
#include <iostream>
#include <map>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "CliqueList.h"

using namespace std;
using namespace boost;

/** Constructor. Note that this isn't particularly fast (but
    doesn't have to be). */
CliqueList::CliqueList(int n, int r, int maxCliqueSize) {
    n_ = n;
    r_ = r;
    max_clique_size_ = maxCliqueSize;
    // number the (hyper)edges
    initEdgeNumbering();
    // initialize the cliques
    clique_.resize(maxCliqueSize+1);
    for(int k=r; k<=maxCliqueSize; k++)
        addCliques(k);
}

/** Adds counts of which cliques are present in a graph.
    Should be multithreaded.
    Possibly not used, as the bookkeeping seems easier if I
    just track the largest clique. (Clearly, if there's a
    9-vertex clique, there are many 8-cliques, and even
    more 7-vertex cliques, etc.) */
void CliqueList::addToCounts(vector<int> count, const bits & g) {
    // loop through the possible sizes of clique
    for(int k = 0; k < int(clique_.size()); k++)
        // loop through the cliques
        for(int i = 0; i < int(clique_[k].size()); i++)
            if (clique_[k][i].is_subset_of(g))
                count[k]++;
}

/** Gets size of the largest clique in a graph.
    Should be multithreaded.
    g: vertices of a hypergraph
    Returns: number of vertices in the largest clique */
int CliqueList::getMaxCliqueSize(const bits & g) {
    // loop through possible clique sizes (starting with largest)
    for(int k = n_; k >= r_; k--)
        // loop through the cliques
        for(int i = 0; i < int(clique_[k].size()); i++) {
            if (clique_[k][i].is_subset_of(g))
                return k;
        }
    return 0;
}

/** Gets vertices covered by hypergraphs.
    g: vertices of a hypergraph
    k: size of the hyperclique to check for
    coveredEdges: ref. to bits, which will be set with
        the edges which are covered
    Side effects: fills in coveredEdges
    Returns: number of k-vertex cliques found */
int CliqueList::getCoveredEdges(const bits & g,
        int k, bits & coveredEdges) {
    // set coveredEdges to all 0 (In theory, if we enforced
    // that callers looked for larger cliques first, we wouldn't
    // need this. But this seems simpler).
    coveredEdges.reset();
    int cliqueCount = 0;
    // loop through the cliques
    for(int i = 0; i < int(clique_[k].size()); i++) {
        if (clique_[k][i].is_subset_of(g)) {
            // record the edges covered by this clique
            coveredEdges |= clique_[k][i];
            ++cliqueCount;
        }
    }
    return cliqueCount;
}

/** Initialize the edge numbering. */
void CliqueList::initEdgeNumbering() {
    int e = 0;
    // this loops through all r-element subsets of n vertices
    SubsetIterator edge(n_, r_);
    do {
        edgeIndex_[ edge.getSet() ] = e;
        ++e;
    } while (edge.next());
}

/** Initialize the cliques with k vertices. */
void CliqueList::addCliques(int k) {
    // this will contain one bitvector per clique
    vector<bits> bitvectors;
    // this loops through the cliques
    SubsetIterator clique(n_, k);
    do {
        // this is the k indices of vertices
        vector<int> &v = clique.getSet();
        // this loops through the edges in a clique
        SubsetIterator edge(k, r_);
        // this tracks the bits corresponding to those edges
        bits edgeBits(edgeIndex_.size());
        do {
            // this will be the r indices of an edge
            // in the clique, among all n vertices
            vector<int> r1;
            vector<int> &e = edge.getSet();
            for(int i=0; i<int(e.size()); i++)
                r1.push_back(v[e[i]]);
            // add this edge into the clique
            edgeBits |= getBitmask(r1);
        } while (edge.next());
        bitvectors.push_back(edgeBits);
    } while (clique.next());
    clique_[k] = bitvectors;                
}

/** Gets the bitmask corresponding to a particular clique. */
bits CliqueList::getBitmask(vector<int> & edge) {
    bits x(edgeIndex_.size());
    // set bits corresponding to each edge
    for(int i=0; i<int(edge.size()); i++)
        x[ edgeIndex_[edge] ] = 1;
    return x;
}


/** Dumps list of cliques to standard output. */
void CliqueList::printCliques() {
    cout << "dumping cliques" << endl;
    for(int cliqueSize = n_; cliqueSize >= r_; --cliqueSize) {
        cout << "n = " << cliqueSize << endl;
        vector<bits> & clique1 = clique_[cliqueSize];
        // loop through cliques of this size
        for(int c = 0; c < int(clique1.size()); ++c) {
            cout << clique1[c] << " ";
        }
        cout << endl;
    }
}

/** Prints the edges in one graph. */
void CliqueList::printGraph(bits &e) {
    // FIXME print the vertices
    cout << e << " ";
}

