/** Various utilities. */

#include <cstdlib>
#include <iostream>
#include <vector>

#include <boost/dynamic_bitset.hpp>

using namespace std;
using namespace boost;

/** Type of bit vector implementation. */
typedef dynamic_bitset<> bits;

/** Creates a random bitset of some size.
    FIXME use a better random number generator */
bits randomBitset(int numBits) {
    bits b(0);
    // first, get enough 64-bit blocks
    for(int i = 0; i < (numBits/64)+1; i++) {
        // get a 64-bit random number
        uint64_t r = random() ^ (random() << 32);
        b.append(r);
    }
    // then truncate to the required size
    b.resize(numBits);
    cout << "made random vector" << endl;
    return b;
}

