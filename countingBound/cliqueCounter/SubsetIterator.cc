
#include "SubsetIterator.h"

#include <vector>

using namespace std;

/** Constructor.
    n: size of set to iterate over
    k: number of things to pick */
SubsetIterator::SubsetIterator(int n, int k) {
    n_ = n;
    k_ = k;
    // initialize to the "first" set
    for(int i=0; i<k_; i++)
        a_.push_back(i);
}

/** Gets the current set. */
vector<int> & SubsetIterator::getSet() {
    return a_;
}

/** Advances to the next set.
    Returns: true iff the next set is defined. If not, it returns
    false (and the value of getSet() is not defined). */
bool SubsetIterator::next() {
    // remove elements of a_ which are "maxed out"
    int max = n_ - 1;
    for(int i = k_-1; i >= 0; --i) {
        if (a_[i] == max)
            a_.pop_back();
        --max;
    }
    // if a_ is empty, everything was "maxed out", and we're done
    if (a_.empty())
        return false;
    // increment the top element of a_
    ++a_.back(); 
    // initialize anything after that
    while (int(a_.size()) < k_)
        a_.push_back( a_.back() + 1 );        
    return true;
}

