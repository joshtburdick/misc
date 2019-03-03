# Optimizing a bound by brute force.

# Given a count of functions, gives a lower bound for at least one
# of the functions (measured in gates). 
# Works with vectors of args.
# Args:
#   m - number of inputs
#   w - number of input wires (in other words, log2(# of functions))
# Returns: minimum number of gates required to implement at least one
#   of the functions.
gates.needed = function(m, w) {
  b = m - 0.5
  ceiling( sqrt( 2*w + b^2 ) - b )
}

# Computes log2 of {n \choose k}, using log scale throughout
# to avoid overflow.
# Args: n, k
# Returns: log2{n \choose k}
choose.log2 = function(n, k) {
  if (length(n) == 1) n = rep(n, length(k))
  if (length(k) == 1) k = rep(k, length(n))
#  stopifnot( length(n) == length(k) )
  stopifnot( all( n >= k ) )
  r = rep(NA, length(n) )
  for(i in 1:length(n)) {
    r[i] = sum(log2( c(((n[i]-k[i]+1):n[i])) )) -
      sum(log2( c(1:k[i]) ))
  }
  r
}

# Lower bound on detecting cliques.
# Args:
#   n, k - number of vertices in the input graph, and the k-cliques
#   s - number of cliques to miss
# Returns: a lower bound on the number of gates.
#   (This is a slightly weaker bound, but should avoid overflow).
bound = function(n, k, s) {

  # number of input wires
  m = choose(n, 2)

  # number of potential cliques (not on log scale)
  num.cliques = choose(n,k) 

  # (log2 of) number of buggy clique-finding functions
  # (This is actually an underestimate; we could miss fewer
  # than s cliques. Presumably, asymptotically, this won't
  # matter much).
  w = choose.log2(num.cliques, s)

  r = gates.needed(m, w) - (s + 1) - 1

  # for now, not clipping this at zero
  r
}

# Naive version of num.functions (more likely to overflow).
num.functions.naive = function(n, k, s)
  log2(choose(choose(n,k), s))

# The clique-detecting bound (naive version, which can and
# does overflow).
# Args:
#   n, k - problem size
#   s.max - different number of cliques to miss
# Returns: the bound, for those values of s.
bound.1 = function(n, k, s) {
  # number of input wires
  m = choose(n, 2)

  # number of cliques we fail to find
  s = c(1:s.max)    # FIXME use e.g. powers of two here?

  # (log2 of) number of buggy circuits
  w = log2( choose(choose(n,k), s) )

  r = gates.needed(m, w) - (s + 1) - 1

  # for now, not clipping this at zero
  r
}


