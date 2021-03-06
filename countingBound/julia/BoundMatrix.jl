#!/home/jburdick/bin/julia
# Writes out the coefficients of the bound matrix, both
# using exact binomials, and using various approximations
# (and log-transformed numbers), as a check of the math.

push!(LOAD_PATH, ".")

using CountingBound
using LogCountingBound

"""
  Writes counts for some values of k, r, and n.
  k: the value of k
  maxN: the maximum value of n
  outputDir: directory in which to write output files
  Side effects: writes files A and b.
    (The choices of r are a bit arbitrary).
"""
function writeCounts(k, maxN, outputDir)
  # values of n to use, from k to maxN
  # (for now, spaced at powers of 2)
  nList = [BigInt(2^i) for i in 1:trunc(log(2,maxN))]
  print(nList)
  print("\n")
  print(typeof(nList[1]))
  print("\n")
  # ??? for some reason, this was giving a syntax error
  # nList = [n for n in nList if n>=k && n<=maxN]
  # therefore, using filter()...
  nList = filter(n -> n>=k && n<=maxN, nList)

  # FIXME should create output directory
  # first, write the bound b (and its log)
  of = open(outputDir * "/b_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,n,bound,exp(log(bound)),log(bound)\n")
  for n in nList
    bound = countingBound(binomial(n, 2), binomial(n, k))
		logBound = NaN
		if bound > 0
			logBound = log(bound)
		end
    write(of, string(k) * "," * string(n) * "," * string(bound) * ","
			* string(exp(logBound)) * ","
      * string(logBound) * "\n")
  end
  close(of)

  # then, write the coefficients
  of = open(outputDir * "/A_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,r,n,A,exp(logA),logA\n")
  for n = nList
    # FIXME what should the bound on r be?
    for r = k:min(n, 2*k)
      print("k=" * string(k) * " r=" * string(r) * " n=" * string(n) * "\n")
      A = approxNumMaximalCliques1(k, r, n)
      logA = logNumMaximalCliques(k, r, n)
      write(of, string(k) * "," * string(r) * "," * string(n)
        * "," * string(A)
				* "," * string(exp(logA))
        * "," * string(logA) * "\n")
    end
  end
  close(of)

end

mkpath("coef1")
# smoke tests
writeCounts(2, 200, "coef1")
writeCounts(3, 200, "coef1")
# almost a smoke test
# writeCounts(4, 100, "coef1")
writeCounts(8, 1000, "coef1")
writeCounts(10, 2000, "coef1")
writeCounts(12, 4000, "coef1")

