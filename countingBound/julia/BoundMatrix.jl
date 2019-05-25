#!/usr/bin/julia
# writes out the coefficients of the bound matrix.

push!(LOAD_PATH, ".")

using CountingBound

# print(countingBound(3,9))
# smoke test
print(approxNumMaximalCliques1(3,4,5))


