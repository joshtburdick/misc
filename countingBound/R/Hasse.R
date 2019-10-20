# Plots a Hasse diagram of selected
# sets of cliques.
# Plots triangular "edges" of a complete graph.

# Three-element subsets of six elements.
vertex.6 = t(combn(6,3))

# Plots a subset of triangles.
#   center: center of where to draw the triangles
#   i: row indices to include
#   r: radius of circle on which to put the points
#   Side effects: draws a subset of triangles
plot.tri.subset = function(center, i, r=1) {
	n = 6
  circle.points = r * cbind(cos(2*pi*(1:n)/n), sin(2*pi*(1:n)/n))
  draw.triangles(circle.points + center,
		vertex.6[i,,drop=FALSE],
    tri.colors$h[i], tri.colors$v[i])
}

pdf("Hasse.pdf")

plot(0,0, c(0,5), c(0,5),
	type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")

plot.tri.subset(c(1,1), c(1), r=0.25)
plot.tri.subset(c(3,3), c(1,2,3), r=0.25)

dev.off()

