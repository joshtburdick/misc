#!/home/jburdick/bin/julia
# writes out the coefficients of the bound matrix.

push!(LOAD_PATH, ".")

using CountingBound

# smoke tests
writeCounts(2, 200, "coef1")
writeCounts(3, 200, "coef1")
# almost a smoke test
# writeCounts(4, 100, "coef1")
# writeCounts(6, 5000000, "coef1")
# writeCounts(8, 1000, "coef1")
# writeCounts(10, 2000, "coef1")
# writeCounts(12, 4000, "coef1")



