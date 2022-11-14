#pragma once

#include <cstdlib>
#include <vector>

#include <boost/dynamic_bitset.hpp>

using namespace std;
using namespace boost;

/** Type of bit vector implementation. */
typedef dynamic_bitset<> bits;

/** Various utilities. */

/** Gets a random bitset. */
bits randomBitset(int numBits);

