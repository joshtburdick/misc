# Functions related to the counting bound. This is similar
# to the original version, but focusses on the case of finding
# k-cliques in graphs with 2k vertices.
# It also uses several approximations:
# - various approximations of binomial coefficients
# - approximates (1-a)^b
# The bound, again, seems to asymptote at a small positive number.

module CountingBound2

export countingBound, approxNumMaximalCliques1, writeCounts

"""
Approximation of binomaial(n, k), when n >> k.
From Wikipedia.
"""
function approxNChooseK(n, k)
  (n/k - 0.5)^k * exp(k) / sqrt(2 * pi * k)
end

"""
Approximation of binomial(2n, n), for n large.
From Wikipedia. (Possibly not used).
"""
function approx2KChooseK(n)
	2 ^ (2*n) / sqrt(pi * n)
end

"""
Counting bound, based on Shannon's argument.
  m: number of edges
  w: number of 'wires' -- that is, log2(number of functions),
    which is the number of bits needed to specify a function
  Returns: average number of NAND gates (with unbounded fan-in)
    required to compute any of those functions.
    (This may not be an integer).
"""
function countingBound(m, w)
  m = BigFloat(m)
  w = BigFloat(w)
  b = m - 0.5
  # the "-1" here is because this is the average, not the max.
  sqrt(2*w + b*b) - b - 1
end

"""
	Bound for finding k-cliques in n-vertex graphs, based on
covering a v-vertex graph with k-hypercliques.
	k: the size of the cliques to find
	n: number of vertices on which we're finding cliques
	v: the size of the larger graph
	Returns: minimal number of NAND gates need to find k-cliques
in an n-vertex graph (based on covering v-vertex graphs
with hypercliques).
"""
function cliqueNandBound(k, n, v)
	# force use of higher-precision arithmetic
	k = BigInt(k)
	n = BigInt(n)
	v = BigInt(v)

	# number of gates needed (for the larger problem)
	numGates1 = countingBound(binomial(v, 2), approxNChooseK(v, k))

	# (log of) probability a k-edge is not covered by an n-clique
	logP = - approxNChooseK(v-k, n-k) / (2^(binomial(n, k)-1))

	# number of k-edges not covered
	numNotCovered = exp(logP) * approxNChooseK(v, k)

	# number of expected n-vertex k-regular hypercliques
	numCliques = approxNChooseK(v, n) / (2^binomial(n, k))

	# to get the bound, we consider the number of gates covered by
	# n-cliques, and then divide by how many of those are expected
	numGates = (numGates1 - numNotCovered) / numCliques

end

"""
  Writes counts for some values of k, r, and n
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
  write(of, "k,n,bound,logBound\n")
  for n in nList
    bound = countingBound(binomial(k, 2), binomial(n, k))
    write(of, string(k) * "," * string(n) * "," * string(bound) * ","
      * string(log(bound)) * "\n")
  end
  close(of)

  # then, write the coefficients
  of = open(outputDir * "/A_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,r,n,A,logA\n")
  for n = nList
    # FIXME what should the bound on r be?
    for r = k:min(n, 2*k)
      print("k=" * string(k) * " r=" * string(r) * " n=" * string(n) * "\n")
      bound = countingBound(binomial(k, 2), binomial(n, k))
      A = approxNumMaximalCliques1(k, r, n)
      logA = logApproxNumMaximalCliques1(k, r, n)
      write(of, string(k) * "," * string(r) * "," * string(n)
        * "," * string(A)
        * "," * string(logA) * "\n")
    end
  end
  close(of)

end

end

