#pragma once

#include <cstdlib>
#include <vector>

#include <boost/dynamic_bitset.hpp>

using namespace std;
using namespace boost;

/** Type of bit vector implementation. */
typedef dynamic_bitset<> bits;

/** Creates a random bitset of some size.
    FIXME
    - use a better random number generator?
    - return this by reference? */
bits randomBitset(int numBits) {
    bits b;
    // first, get enough 64-bit blocks
    for(int i = 0; i < (numBits/64)+1; i++) {
        // get a 64-bit random number
        uint64_t r = random() ^ (random() << 32);
        b.append(r);
    }
    // then truncate to the required size
    b.resize(numBits);
    return b;
}
`

/** Various utilities. */

/** Gets a random bitset. */
bits randomBitset(int numBits);

