# Simple tests of some bounding functions.

source("gateCount.R")

# log2(choose)
{
  set.seed(1)
  k = sample(20, 100, replace=T) 
  n = k + sample(20, 100, replace=T)
  c1 = choose(n, k)
  c2 = 2 ^ choose.log2(n, k)
  cat("max. abs. choose(n,k) diff = ", max(abs(c1-c2)), "\n")
}



