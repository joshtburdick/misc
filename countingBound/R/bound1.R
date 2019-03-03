# Bound based on probability of seeing a clique of a given size.
# Most of these functions accept either a single number,
# or equal-length vector arguments.

# Counting bound, based on Shannon's argument.
#   m: number of edges
#   w: number of "wires" -- that is, log2(number of functions),
#     which is the number of bits needed to specify a function
#   Returns: average number of gates required to compute
#     any of those functions. (This may not be an integer).
counting.bound = function(m, w) {
  b = m - 0.5
  # the "-1" here is because this is the average, not the max.
  sqrt(2*w + b*b) - b - 1
}

# Computes a bound for finding k-cliques.
# Note that k <= r <= n.
#   k: size of cliques sought
#   r: number of vertices in "clique detector" whose minimal
#     size we're inferring
#   n: number of vertices in the larger graph (which can
#     be a vector)
# Returns: number of gates for each value of n
approx.bound.1 = function(k, r, n) {
  # number of gates (from the Shannon bound)
  min.total.gates = counting.bound(choose(k, 2), choose(n, k))

  # expected number of edges "left over" (not covered
  # by an r-clique)
  p = (1 - 2 ^ (-(choose(r, k)-1))) ^ choose(n-k, r-k)
  num.left.over = (p/2) * choose(n, k)

  # expected number of r-cliques
  num.r.cliques = choose(n, r) * (2 ^ (-choose(r, k)))

  # the remaining gates are part of some number of
  # clique detectors
  r = (min.total.gates - num.left.over) / num.r.cliques

  r
}

#### versions of these using logs
# (deprecated; switching to using Julia exact arithmetic instead)

# Like approx.bound.1, but using logs to avoid overflow.
# (Hopefully this is using logs in enough places; for now,
# not doing quite so many things in log-space).
approx.bound.2 = function(k, r, n) {
  # number of gates (from the Shannon bound)
  # not log-transformed (may or may not be a problem; hopefully
  # not, as there's no 2^___ here)
  min.total.gates = counting.bound(choose(k, 2), choose(n, k))

  # expected number of edges "left over" (not covered
  # by an r-clique)
#  was: p = (1 - 2 ^ (-(choose(r, k)-1))) ^ choose(n-k, r-k)
# note that (I think) 2^blah = exp(blah / log(2))
  # also using the bound based on
  #   lim a->\inf (1 - 1/a)^a = e^-1
#  (hopefully "choose(r,k)-1" won't be too large)
  l.p = - exp( lchoose(n-k, r-k) - log(2)*log(choose(r,k)-1) )
  l.num.left.over = l.p + lchoose(n, k) - log(2)

  # expected number of r-cliques
  # was: num.r.cliques = choose(n, r) * (2 ^ (-choose(r, k)))
  l.num.r.cliques = lchoose(n, r) - (log(2)*lchoose(r, k))

  # the remaining gates are part of some number of
  # clique detectors
#  l.r = log(min.total.gates - exp(l.num.left.over)) - l.num.r.cliques
#  exp(l.r)

  l.r = log(min.total.gates - exp(l.num.left.over)) - l.num.r.cliques
  exp(l.r)
}

# Computes an approximate lower bound on finding k-cliques
# in r-vertex graphs.
#   k: size of the clique we're looking for
#   r: size of the larger circuit (which can be a vector)
#   Returns: smallest number of NAND gates needed to find
#     k-cliques in r-vertex graphs.
approx.bound = function(k, r) {
  n = c(r:(3*r))
  g = approx.bound.1(k, r, n)
  cbind(n, g)
}

# z = approx.bound(6, 10)

