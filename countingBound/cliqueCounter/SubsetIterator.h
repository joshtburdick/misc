#pragma once

#include <vector>

using namespace std;

/** Iterates through subsets of a set of numbers
    (in lexicographic order). */
class SubsetIterator {

public:
    /** Constructor.
        n: size of set to iterate over
        k: number of things to pick */
    SubsetIterator(int n, int k);

    /** Gets the current set. */
    vector<int> & getSet();

    /** Advances to the next set.
        Returns: true iff the next set is defined. If not, it returns
        false (and the value of getSet() is not defined). */
    bool next();

private:
    int n_;

    int k_;

    vector<int> a_;
};

