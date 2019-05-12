#!/usr/bin/julia
# Looks for a bound.
# I'm not sure if it's currently working.

"""
Ternary search (based on Wikipedia). This searches for an integer
    value maximizing f.
  f: function to maximize
  left, right: bounds
  Returns: value maximizing f
"""
function ternarySearch(f, left, right)
  print(string(left) * " " * string(right) * "\n")
  # left and right are the current bounds; 
  # the maximum is between them
  if abs(right - left) <= 1
    return trunc((left+right)/2)
  end
  leftThird = ((2*left + right)/3)
  rightThird = ((left + 2*right)/3)
  if f(trunc(leftThird)) < f(trunc(rightThird))
    return ternarySearch(f, leftThird, right) 
  else
    return ternarySearch(f, left, rightThird)
  end
end

"""
Counting bound, based on Shannon's argument.
  m: number of edges
  w: number of 'wires' -- that is, log2(number of functions),
    which is the number of bits needed to specify a function
  Returns: average number of gates required to computea
    any of those functions. (This may not be an integer).
"""
function countingBound(m, w)
  b = m - 0.5
  # the "-1" here is because this is the average, not the max.
  sqrt(2*w + b*b) - b - 1
end

"""
  Number of maximal hypercliques of some size.
  Note that k < r < n .
  k: number of vertices per hyperedge
  r: number of vertices in the clique
  n: vertices in the larger graph
  Returns: number of maximal hypercliques
"""
function numMaximalCliques(k, r, n)
  k = BigInt(k)
  r = BigInt(r)
  n = BigInt(n)
  one = BigInt(1)
  two = BigInt(2)

  # expected number of r-cliques
  numRCliques = Rational(binomial(n, r),
    two ^ binomial(r, k))

  # probability that one of those is not covered by a larger clique
  a = two ^ binomial(r, k-one)
  p = Rational(a-one, a) ^ (n-r)
  # println(a)
  # println(p)
  # result is number of cliques, * prob. they're maximal
  numRCliques * p
end

"""
  Bound, using one size of clique.
    This is approximate, because it's approximating the number of
    edges and cliques.
  k: number of vertices per hyperedge
  r: number of vertices for the clique finder
  n: vertices in the larger graph
  Returns: number of NAND gates
"""
function bound2(k, r, n)
  k = BigInt(k)
  r = BigInt(r)
  n = BigInt(n)
  one = BigInt(1)
  two = BigInt(2)

  # the bound on the total number of gates
  minTotalGates = ceil(countingBound(binomial(k, two), binomial(n, k)))

  # expected number of edges "left over" (not covered
  # by an r-clique)
  # p = (1 - 2 ^ (-(choose(r, k)-1))) ^ choose(n-k, r-k)
  a = two ^ (binomial(r, k) - one)
  # p = Rational(a-one, a) ^ binomial(n-k, r-k)
  # print("0\n")
  b = binomial(n-k, r-k)
  p = Rational((a-one)^b, a^b)
  # print("1\n")
  # num.left.over = (p/2) * choose(n, k)
  numLeftOver = p * binomial(n, k) / two
  # print("2\n")

  # expected number of r-cliques
  numRCliques = Rational(binomial(BigInt(n), BigInt(r)),
    BigInt(2) ^ binomial(r, k))
  # print("3\n")
  # the remaining gates are part of some number of
  # clique detectors
  r = (minTotalGates - numLeftOver) / numRCliques
end

"""
  Bound, using one size of clique.
    This is approximate, as it's estimating the number of edges
    and cliques. It also approximates the "(1 - 1/a)^..." bit.
  k: number of vertices per hyperedge
  r: number of vertices for the clique finder
  n: vertices in the larger graph
  Returns: number of NAND gates
"""
function bound3(k, r, n)
  k = BigInt(k)
  r = BigInt(r)
  n = BigInt(n)
  one = BigInt(1)
  two = BigInt(2)

  # the bound on the total number of gates (in the larger graph)
  minTotalGates = ceil(countingBound(binomial(k, two), binomial(n, k)))

  # expected number of edges "left over" (not covered by an r-clique)
  # First, the probability of a hyperedge being "left over"
  # (not covered by an r-clique) is
  #
  # p = (1 - 2 ^ (-(choose(r, k)-1))) ^ choose(n-k, r-k)
  #
  # We use an approximate upper bound (which should be good
  # for large numbers).
  logP = - binomial(n-k, r-k) / (two ^ (binomial(r, k) - one))
  # (log of) number of edges left over
  logNumLeftOver = logP + log(binomial(n, k) / two)

  # expected number of r-cliques
  # (Using right-shift for division, and adding 1 in case
  # of rounding error)
  numRCliques = (binomial(n, r) >> binomial(r, k)) + 1
  # the remaining gates are part of some number of
  # clique detectors
  r = (minTotalGates - exp(logNumLeftOver)) / numRCliques
end

function printBound1(k, r, n)
  b = bound1(k,r,n)
  print(string(k, ",", r, ",", n, ",", b, "\n"))
end

# print("k,r,n,bound\n")
# for n in 100:20:600
#   printBound1(5, 6, n)
# end




