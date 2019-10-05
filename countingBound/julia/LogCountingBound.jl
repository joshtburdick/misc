# Functions related to the counting bound.
# This version uses analytic approximations from Wikipedia,
# and works with log-transformed numbers.
# (All logs are natural, unless otherwise noted).

module LogCountingBound

export logCountingBound, logApproxNumMaximalCliques,
  approxLogNChooseK, approxBoundAlmostOneExp, logDiff
	# ??? export fewer of these?

"""
Approximate log of binomial(n, k), when n >> k.
Based on Wikipedia.
FIXME use better approximation from Wikipedia?
"""
function approxLogNChooseK(n, k)
  k * log(n/k - 0.5) + k - 0.5 * log(2 * pi * k)
end

"""
Upper-bound on

log( (1 - 1/a)^b )

given the log of a and b. This is more accurate when
a is large, and b is small.
  logA, logB: log of a and b, respectively
  Returns: log of the above expression
"""
function approxBoundAlmostOneExp(logA, logB)
  - exp( logB - logA )
end

"""
Given log(a) and b, computes log( a +/- b ).
	logA: a log-transformed number
	b: a (non-log-transformed) number, smaller than e^a
	plus: if true, add b; if false, subtract b
	Returns: log(a +/- b), or NaN if something goes negative.
(This is computed using techniques adapted from those used
to deal with log-likelihoods e.g. in NLP).
"""
function logPlusMinus(logA, b, plus)
	if b <= 0
		return NaN
	end
	# compute b, as a fraction of a
	bScaled = exp( log(b) - logA )
	# add or subtract that from a
	diff = if plus; 1 + bScaled; else 1 - bScaled; end
	if diff <= 0
		return NaN
	end
	# return log of that, back on original scale
	log( diff ) + logA
end

"""
Counting bound, based on Shannon's argument.
  m: number of edges
  logW: log of number of 'wires' -- that is,
		log(log2(number of functions)), where
		'log2(number of functions)' is the number of bits
		needed to specify a function
  Returns: log of average number of NAND gates (with unbounded
		fan-in) required to compute any of those functions.
    (This may not be an integer).
"""
function logCountingBound(m, logW)
  m = BigFloat(m)
  logW = BigFloat(logW)
  b = m - 0.5
	# this presumably should be this:
  # the "-1" here is because this is the average, not the max.
  # bound = sqrt(2*w + b*b) - (b - 1)
	# but we use log-transformed numbers in several places
	x = logPlusMinus(log(2) + logW, b*b, true) / 2
	bound = logPlusMinus(x, b - 1, false)
	bound
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
  logA = log(BigInt(2)) * binomial(r, k-1)
	logP = approxBoundAlmostOneExp(logA, log(n-r))

	logNumMaximalCliques = (
    # number of possible cliques of this size
		approxLogNChooseK(n, r)
    # P(clique is present ...)
		+ (- log(BigInt(2)) * binomial(r, k))
    # P( ... and not covered by a larger clique)
    + logP)

  logNumMaximalCliques
end

end

