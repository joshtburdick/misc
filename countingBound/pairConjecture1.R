
# Args:
#   n - number of pairs of numbers
# Returns: matrix with two rows, 
random.pairs = function(n) {
  a = matrix(sample(2*n), nrow=2)
  apply(a, 2, sort)
}

# Plots the distribution of the difference between the
# larger and smaller number in each pair.
# Args:
#   n - number of pairs of numbers
# Side effects: plots a histogram
plot.random.pair.diff = function(n) {
  a = random.pairs(n)
  hist(a[2,] - a[1,], breaks=100,
    col="lightgrey", lwd=0.1,
    main=paste("n =", n))
  abline(v = n, col="#0000ffa0")
}

set.seed(0)
pdf("pairConjecture1.pdf", width=7.5, height=7.5)
par(mfrow=c(3,3))
for(n in c(10, 50, 100, 1e3, 1e3, 1e3, 1e4, 1e4, 1e5))
  plot.random.pair.diff(n)
dev.off()

