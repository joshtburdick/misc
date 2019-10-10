#!/home/jburdick/bin/julia
# Writes out the coefficients of the bound matrix, both
# using exact binomials, and using various approximations
# (and log-transformed numbers), as a check of the math.

push!(LOAD_PATH, ".")

using LogCountingBound

"""
  Writes counts for some values of k, r, and n.
  k: the value of k
  maxN: the maximum value of n (this will use
		powers of two	<= maxN)
	logADropOff: this will start r at k, and continue as the log bound
		increases, until the bound drops below (max - logADropOff)
  outputDir: directory in which to write output files
  Side effects: writes files A and b.
    (The choices of r are a bit arbitrary).
"""
function writeCounts(k, maxN, logADropOff, outputDir)
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
		# tracks the maximum (initially undefined)
		maxLogA = NaN
		# loop through possible values of r (this will usually stop early)
    for r = k:n
      print("k=" * string(k) * " r=" * string(r) * " n=" * string(n) * "\n")
			# compute approximation, and print it
      logA = logApproxNumMaximalCliques(k, r, n)
      write(of, string(k) * "," * string(r) * "," * string(n)
        * "," * string(logA) * "\n")
			# keep track of maximum
			if isnan(maxLogA) || logA > maxLogA
				maxLogA = logA
			end
			# if we've fallen far enough below the maximum, stop the loop
			if logA <= maxLogA - logADropOff
				break
			end
			# ??? sometimes this seems to stop after one number,
			# when r is very large
    end
  end
  close(of)

end

logADropOff = 100 * log(10)

outputDir = "logCoef1"
mkpath(outputDir)

writeCounts(6, 1e100, logADropOff, outputDir)
writeCounts(8, 1e100, logADropOff, outputDir)
writeCounts(10, 1e100, logADropOff, outputDir)
writeCounts(12, 1e100, logADropOff, outputDir)
writeCounts(20, 1e100, logADropOff, outputDir)
writeCounts(30, 1e100, logADropOff, outputDir)
writeCounts(100, 1e100, logADropOff, outputDir)

writeCounts(1000, 1e100, logADropOff, outputDir)

writeCounts(10000, 1e100, logADropOff, outputDir)

