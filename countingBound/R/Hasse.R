# Plots a Hasse diagram of selected
# sets of cliques.
# Plots triangular "edges" of a complete graph.

# Three-element subsets of six elements.
vertex.6 = t(combn(6,3))

x1 = dnorm(0:20, mean=10, sd=5)
x1 = x1 / max(x1)

# Transforms from 3-D to screen coordinates.
# x: the user coordinates, as a three-row matrix
#		(or a column vector), with rows:
#   x: in [-1,1]
#     FIXME: scale to Hasse diagram width?
#   y: in [0,20], as it's the number of triangles present
#   z: in [0,10]: this is arbitrarily heights (for illustration;
#     they're presumably inaccurate)
p = function(x) {
	A = rbind(c(7, 0.3, 0), c(0, 0.2, 1))
#	A = rbind(c(1, 0, 0), c(0,0,1))
	A %*% x
}

# Plots counts, sort of.
plot.counts = function() {
	b = rbind(x1, 0:20, 0)	
	p1 = t(p(b))
	lines(p1)

	b = rbind(-x1, 0:20, 0)	
	p1 = t(p(b))
	lines(p1)
}

# Plots a subset of triangles.
#   center: center of where to draw the triangles
#   i: row indices to include
#		label: label for the set
#   r: radius of circle on which to put the points
#   Side effects: draws a subset of triangles
plot.tri.subset = function(center, i, label, r=1) {
	n = 6
  circle.points = r * cbind(cos(2*pi*(1:n)/n), sin(2*pi*(1:n)/n))
#  print(circle.points)
	p = t(t(circle.points) + center)
  draw.triangles(p,
		vertex.6[i,,drop=FALSE],
    tri.colors$h[i], tri.colors$v[i])
}

# Wrapper for the above which:
#   - uses the projection transformation
#   - draws a line down to the "ground"
#   x, z: coordinates for the set
#   i: the set (which will determine the y coordinate)
plot.tri.subset.3d = function(x, z, i) {
	# create the coordinates
	y = length(unique(i))
	p0 = as.vector(p(t(t(c(x, y, 0)))))
	p1 = as.vector(p(t(t(c(x, y, z)))))
	lines(rbind(p0, p1), col="#00000040", lwd=3)
	# for cheap perspective, make distant sets a bit smaller
	plot.tri.subset(p1, i, r=1)  # -z/50)
}

pdf("Hasse.pdf")
# par(mar=c(0,0,0,0))
world.bounds = cbind(c(-1,0,0), c(1,20,10))
screen.bounds = range(p(world.bounds))

plot(0,0,
#	xlim=range(screen.bounds[1,]), ylim=range(screen.bounds[2,]),
	xlim = screen.bounds, ylim = screen.bounds,
	type="n")   #, xaxt="n", yaxt="n",xlab="", ylab="", bty="n")

plot.counts()

plot.tri.subset.3d(0, 1, c(1))
# plot.tri.subset.3d(0, 3, c(1,2))
plot.tri.subset.3d(0, 3, c(1,2,5,11))
v = which(apply(vertex.6, 1, function(a) all(a!=6)))
plot.tri.subset.3d(0, 5, v)
plot.tri.subset.3d(0, 7, c(1:20))
dev.off()

