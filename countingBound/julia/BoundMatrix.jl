#!/usr/bin/julia
# writes out the coefficients of the bound matrix.

push!(LOAD_PATH, ".")

using CountingBound


# smoke test
print(approxNumMaximalCliques1(6,12,10000000))
print("\n")

