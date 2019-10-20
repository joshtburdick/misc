#!/home/jburdick/bin/julia
# Writes out the coefficients of the bound matrix, both
# using exact binomials, and using various approximations
# (and log-transformed numbers), as a check of the math.

push!(LOAD_PATH, ".")

using CountingBound, LogCountingBound

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
  # print(nList)
  # print("\n")
  # print(typeof(nList[1]))
  # print("\n")
  # ??? for some reason, this was giving a syntax error
  # nList = [n for n in nList if n>=k && n<=maxN]
  # therefore, using filter()...
  nList = filter(n -> n>=k && n<=maxN, nList)

  # FIXME should create output directory
  # first, write the bound b
  of = open(outputDir * "/b_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,n,logBound\n")
  for n in nList
    bound = countingBound(binomial(n, 2), binomial(n, k))
		logBound = NaN
		if bound > 0
			logBound = log(bound)
		end
    write(of, string(k) * "," * string(n) * ","
      * string(logBound) * "\n")
  end
  close(of)

  # then, write the coefficients
  of = open(outputDir * "/A_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,r,n,logA\n")
  for n = nList
		# tracks the maximum (initially undefined)
		maxLogA = NaN
		# loop through possible values of r (this will usually stop early)
    for r = k:n
      print("k=" * string(k) * " r=" * string(r) * " n=" * string(n) * "\n")
			# compute log of number of cliques, and print it
      logA = logNumMaximalCliques(k, r, n)
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
    end
  end
  close(of)

end

outputDir = "logCoef1"
mkpath(outputDir)

for k in [6,8,10,12,14,16,18,20,50,100,200,500,1000]
	writeCounts(k, 1e60, 100 * log(10), outputDir)
end

