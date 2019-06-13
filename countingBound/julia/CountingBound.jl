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
  m = BigFloat(m)
  w = BigFloat(w)
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
  Log of the approximate number of maximal hypercliques of some size.
  Note that k < r < n .
  Also, the precision of what's returned can be set by setprecision().
  k: number of vertices per hyperedge
  r: number of vertices in the clique
  n: vertices in the larger graph
  Returns: (natural log of) the number of maximal hypercliques. This is
    approximate, because it's the expected number, and it uses an
    approximation for (1 - 1/a)^b. (Again, presumably it's more
    accurate for larger numbers).
"""
function logApproxNumMaximalCliques1(k, r, n)
  k = BigInt(k)
  r = BigInt(r)
  n = BigInt(n)
  one = BigInt(1)
  two = BigInt(2)

  # probability that one of those is not covered by a larger clique
  a = log(2) * binomial(r, k-1)
  logP = a / (n-r)

  logR = (log(binomial(n, r))  # number of possible cliques
    - log(2) * binomial(r, k)  # P(clique is present) ...
    + logp                     # P(and not part of a larger clique)
  logR 
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
  write(of, "k,n,bound\n")
  for n in nList
    bound = countingBound(binomial(k, 2), binomial(n, k))
    write(of, string(k) * "," * string(n) * "," * string(bound) * "\n")
  end
  close(of)

  # then, write the coefficients
  of = open(outputDir * "/A_k=" * string(k) * "_maxN=" * string(maxN) * ".csv", "w")
  write(of, "k,r,n,A\n")
  for n = nList
    # FIXME what should the bound on r be?
    for r = k:min(n, 2*k)
      print("k=" * string(k) * " r=" * string(r) * " n=" * string(n) * "\n")
      bound = countingBound(binomial(k, 2), binomial(n, k))
      A = approxNumMaximalCliques1(k, r, n)
      write(of, string(k) * "," * string(r) * "," * string(n) * "," * string(A) * "\n")
    end
  end
  close(of)

end

end

