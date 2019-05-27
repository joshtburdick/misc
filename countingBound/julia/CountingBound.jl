# Functions related to the counting bound.

module CountingBound

export countingBound, approxNumMaximalCliques1, writeCounts

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
  b = m - 0.5
  # the "-1" here is because this is the average, not the max.
  sqrt(2*w + b*b) - b - 1
end

"""
  Approximate number of maximal hypercliques of some size.
  Note that k < r < n .
  Also, the precision of what's returned can be set by setprecision().
  k: number of vertices per hyperedge
  r: number of vertices in the clique
  n: vertices in the larger graph
  Returns: number of maximal hypercliques. This is approximate, because
    it's the expected number. (Presumably it's more accurate for larger
    numbers).
"""
function approxNumMaximalCliques1(k, r, n)
  k = BigInt(k)
  r = BigInt(r)
  n = BigInt(n)
  one = BigInt(1)
  two = BigInt(2)

  # probability that one of those is not covered by a larger clique
  a = one << binomial(r, k-one)
  print("computed a\n")
  pNumerator = (a-one) ^ (n-r)
  print("computed numerator\n")
  pDenominator = a ^ (n-r)
  print("computed denominator\n")

  # expected number of r-cliques should be equivalent to:
  # numRCliques = Rational(binomial(n, r), two ^ binomial(r, k))
  # # result is number of cliques, * prob. they're maximal
  # numRCliques * (pNumerator / pDenominator)
  r = (pNumerator * binomial(n, r)) /
    (pDenominator * (one << binomial(r, k)))
  # ??? is this off by two?
  r
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
  # FIXME should create output directory
  # first, write the bound b
  of = open(outputDir * "/b_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,n,bound\n")
  for n = k:maxN
    bound = countingBound(binomial(k, 2), binomial(n, k))
    write(of, string(k) * "," * string(n) * "," * string(bound) * "\n")
  end
  close(of)

  # then, write the coefficients
  of = open(outputDir * "/A_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,r,n,A\n")
  for n = k:maxN
    for r = k:min(n, 2*k)
      bound = countingBound(binomial(k, 2), binomial(n, k))
      A = approxNumMaximalCliques1(k, r, n)
      write(of, string(k) * "," * string(r) * "," * string(n) * "," * string(A) * "\n")
    end
  end
  close(of)

end

end

