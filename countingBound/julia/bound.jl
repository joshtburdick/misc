
"""
Counting bound, based on Shannon's argument.
  m: number of edges
  w: number of 'wires' -- that is, log2(number of functions),
    which is the number of bits needed to specify a function
  Returns: average number of gates required to compute
    any of those functions. (This may not be an integer).
"""
function countingBound(m, w)
  b = m - 0.5
  # the "-1" here is because this is the average, not the max.
  sqrt(2*w + b*b) - b - 1
end

"""
  Bound (approximate, because it's approximating the number of
    edges and cliques).
  k: number of vertices per hyperedge
  r: number of vertices for the clique finder
  n: vertices in the larger graph
  Returns: number of NAND gates
"""
function bound1(k, r, n)
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

function printBound1(k, r, n)
  b = bound1(k,r,n)
  print(string(k, ",", r, ",", n, ",", b, "\n"))
end

print("k,r,n,bound\n")
for n in 100:20:600
  printBound1(5, 6, n)
end

