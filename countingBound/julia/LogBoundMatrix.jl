#!/home/jburdick/bin/julia
# Writes out the coefficients of the bound matrix, both
# using exact binomials, and using various approximations
# (and log-transformed numbers), as a check of the math.

push!(LOAD_PATH, ".")

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
  # first, write the bound b
  of = open(outputDir * "/b_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,n,approxLogBound\n")
  for n in nList
		logW = approxLogNChooseK(n, k)
		logBound = logCountingBound(binomial(n, 2), logW)
    write(of, string(k) * "," * string(n) * ","
      * string(logBound) * "\n")
  end
  close(of)

  # then, write the coefficients
  of = open(outputDir * "/A_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,r,n,approxLogA\n")
  for n = nList
    # FIXME what should the bound on r be?
    for r = k:min(n, 2*k)
      print("k=" * string(k) * " r=" * string(r) * " n=" * string(n) * "\n")
      logA = logApproxNumMaximalCliques(k, r, n)
      write(of, string(k) * "," * string(r) * "," * string(n)
        * "," * string(logA) * "\n")
    end
  end
  close(of)

end

outputDir = "logCoef1"
mkpath(outputDir)
# smoke tests
writeCounts(2, 200, outputDir)
writeCounts(3, 200, outputDir)
writeCounts(8, 1000, outputDir)
writeCounts(10, 2000, outputDir)
writeCounts(12, 4000, outputDir)

