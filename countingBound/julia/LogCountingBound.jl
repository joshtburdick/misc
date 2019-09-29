# Functions related to the counting bound.
# This version uses analytic approximations from Wikipedia,
# and works with log-transformed numbers.
# (All logs are natural, unless otherwise noted).

module LogCountingBound

export logCountingBound, logApproxNumMaximalCliques

"""
Approximate log of binomial(n, k), when n >> k.
Based on Wikipedia.
"""
function approxLogNChooseK(n, k)
  k * log(n/k - 0.5) + k - 0.5 * log(2 * pi * k)
end

"""
Upper-bound on

log( (1 - 1/a)^b )

given the log of a and b. This is more accurate when
a and b are large.
  logA, logB: log of a and b, respectively
  Returns: log of the above expression
"""
function approxBoundAlmostOneExp(logA, logB)
  - exp( logB - logA )
end

"""
Counting bound, based on Shannon's argument.
  m: number of edges
  w: number of 'wires' -- that is, log2(number of functions),
    which is the number of bits needed to specify a function
  Returns: average number of NAND gates (with unbounded fan-in)
    required to compute any of those functions.
    (This may not be an integer).
??? compute this log-transformed?
"""
function logCountingBound(m, w)
  m = BigFloat(m)
  w = BigFloat(w)
  b = m - 0.5
  # the "-1" here is because this is the average, not the max.
	# FIXME use an approximation here (in case m is humungous) ?
  bound = sqrt(2*w + b*b) - b - 1
	if bound >= 1
		return(log(bound))
	else
		return(0)
	end
end

"""
  Log of the approximate number of maximal hypercliques of some size.
  Note that k < r < n .
  Also, the precision of what's returned can be set by setprecision().
  k: number of vertices per hyperedge
  r: number of vertices in the clique
  n: vertices in the larger graph
  Returns: (natural log of) the number of maximal hypercliques. This
    is approximate, because it's the expected number, and it uses an
    approximation for (1 - 1/a)^b. (Again, presumably it's more
    accurate for larger numbers).
  We assume that binomial(r, k-1) and binomial(r, k) are
		computable exactly, but that n may be much larger, and so
		'n choose r' or 'n choose k' require approximation.
"""
function logApproxNumMaximalCliques(k, r, n)
  k = BigInt(k)
  r = BigInt(r)
  n = BigInt(n)

  # log of probability that one of those is
	# not covered by a larger clique
  a = log(BigInt(2)) * binomial(r, k-1)
	logP = approxBoundAlmostOneExp(a, log(n-r))

	logNumMaximalCliques =
    # number of possible cliques of this size
		approxLogNChooseK(n, r)
    # P(clique is present ...)
		+ (- log(BigInt(2)) * binomial(r, k))
    # P( ... and not covered by a larger clique)
    + logP 

  logNumMaximalCliques
end

end

