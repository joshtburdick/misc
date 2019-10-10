# Functions related to the counting bound.
# This version uses log-transformed numbers (but doesn't
# use analytic approximations to binomial coefficients).
# All logs are natural, unless otherwise noted.

module LogCountingBound

export logApproxNumMaximalCliques

"""
  Log of the approximate number of maximal hypercliques of some size.
  Note that k < r < n .
  Also, the precision of what's returned can be set by setprecision().
  k: number of vertices per hyperedge
  r: number of vertices in the clique
  n: vertices in the larger graph
  Returns: (natural log of) the number of maximal hypercliques.
		This is approximate, as it's the expected number, and multiplies
		floating point by adding logarithms.
"""
function logNumMaximalCliques(k, r, n)
  k = BigInt(k)
  r = BigInt(r)
  n = BigInt(n)

  # log of probability that one of those is
	# not covered by a larger clique
  a = BigInt(2) ^ binomial(r, k-1)
	# ??? will this have enough significant figures?
	logP = (log(a-1) - log(a)) * (n - r)

	logNumMaximalCliques = (
    # number of possible cliques of this size
		log(binomial(n, r))
    # P(clique is present ...)
		+ (- log(BigInt(2)) * binomial(r, k))
    # P( ... and not covered by a larger clique)
    + logP)

  logNumMaximalCliques
end

end

