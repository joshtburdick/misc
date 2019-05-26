# Functions related to the counting bound.

module CountingBound

export countingBound, approxNumMaximalCliques1

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
    (two * pDenominator * (one << binomial(r, k)))

end



end

