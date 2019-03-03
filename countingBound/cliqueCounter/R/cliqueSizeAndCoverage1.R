# Plots counts of clique size and coverage.

library(ggplot2)
library(reshape2)

# Reads in one file of counts.
read.counts = function(n, k) {
  f = paste0("../counts2/coverage_", n, "_", k, "_7_1000.txt")
  a = read.table(f, as.is=T)
  m = ncol(a) / 2

  num.cliques = a[,1:m]
  colnames(num.cliques) = 0:(m-1)
  num.cliques = cbind(n, k, melt(num.cliques, id.vars=c()))
  colnames(num.cliques) = c("n", "k", "r", "num.cliques")

  edges.covered = a[,(m+1):(2*m)]
  # uncomment this to show relative amounts covered
#  edges.covered = edges.covered / apply(edges.covered, 1, max)
  colnames(edges.covered) = 0:(m-1)
  edges.covered = cbind(n, k, melt(edges.covered, id.vars=c()))
  colnames(edges.covered) = c("n", "k", "r", "edges.covered")

  list(num.cliques=num.cliques, edges.covered=edges.covered)
}

# read in counts for various n and k
counts = list()
for(n in c(7:25))
  for(k in c(2:4))
    counts[[paste0(n,"_",k)]] = read.counts(n,k)
num.cliques = do.call(rbind, lapply(counts, function(x) x$num.cliques))
edges.covered = do.call(rbind, lapply(counts, function(x) x$edges.covered))
edges.covered$edges.covered.frac = edges.covered$edges.covered /
  choose(edges.covered$n, edges.covered$k)

# restrict to cases with num. vertices >= num. vertices in an edge
# XXX the use of "factors" makes this complicated; I'm sure there's
# a cleaner way, but I don't know what it is
num.cliques = num.cliques[
  as.numeric(as.character(num.cliques$r)) >= num.cliques$k , ]
edges.covered = edges.covered[
  as.numeric(as.character(edges.covered$r)) >= edges.covered$k , ]
# shuffle rows of these (to hopefully clarify display)
i = sample(nrow(num.cliques))
num.cliques = num.cliques[i,]
edges.covered = edges.covered[i,]

# Expected number of cliques of size r in a
# k-regular n-vertex hypergraph.
exp.num.cliques = function(n, k, r) {
  r = choose(n, r) * (2 ^ (-choose(r, k)))
  # if the number of cliques is less than 1, set it to 1
  r[ r < 1 ] = 1
  r
}

# Expected number of hyperedges covered by
# cliques of size r in a k-regular n-vertex hypergraph
exp.num.edges.covered = function(n, k, r) {
  p = 1 - (1 - 2 ^ (-(choose(r, k)-1))) ^ choose(n-k, r-k)
  r = (p/2) * choose(n, k)
  r[ r < 1 ] = 1
  r
}

e = unique(num.cliques[,c("n","k","r")])
e$num.cliques = exp.num.cliques(e$n, e$k, as.numeric(as.character(e$r)))
e$edges.covered = exp.num.edges.covered(e$n, e$k, as.numeric(as.character(e$r)))
e$edges.covered.frac = e$edges.covered / choose(e$n, e$k)
expected = e
rm(e)


png("numCliques.png", width=7.5, height=4, units="in", res=150)
f = ggplot(num.cliques, aes(n, num.cliques, color=r)) +
  geom_point(alpha=0.01, size=0.2) +
  geom_smooth(size=0.2, alpha=1, se=FALSE) +
  facet_grid(. ~ k) + scale_y_log10() + theme_bw() +
  geom_smooth(data=expected, se=FALSE, linetype=2, size=0.6, alpha=1)
# FIXME this darkens stuff, but applies to the whole graph
# + scale_color_hue(l=40, c=30)
show(f)
dev.off()

png("edgesCovered.png", width=7.5, height=4, units="in", res=150)
f = ggplot(edges.covered, aes(n, edges.covered, color=r)) +
  geom_point(alpha=0.01, size=0.2) +
  geom_smooth(size=0.2, alpha=1, se=FALSE) +
  facet_grid(. ~ k) + scale_y_log10() + theme_bw() +
  geom_smooth(data=expected, se=FALSE, linetype=2, size=0.6, alpha=1)
show(f)
dev.off()

png("fracCovered.png", width=7.5, height=4, units="in", res=150)
f = ggplot(edges.covered, aes(n, edges.covered.frac, color=r)) +
  geom_point(alpha=0.01, size=0.2) +
  geom_smooth(size=0.2, alpha=1, se=FALSE) +
  facet_grid(. ~ k) + theme_bw() +
  geom_smooth(data=expected, se=FALSE, linetype=2, size=0.6, alpha=1)
show(f)
dev.off()


if (F) {
png("fracCovered.png", width=7.5, height=3, units="in", res=150)
f = ggplot(expected, aes(n, edges.covered.frac, color=r)) +
  geom_smooth() +
  facet_grid(. ~ k) + theme_bw()  # + scale_y_log10()
show(f)
dev.off()
}
