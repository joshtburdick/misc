#!/usr/bin/julia
# writes out the coefficients of the bound matrix.

push!(LOAD_PATH, ".")

using CountingBound

# smoke test
writeCounts(2, 20, "coef1")


