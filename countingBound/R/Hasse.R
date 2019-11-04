# Plots a Hasse diagram of selected
# sets of cliques.
# Plots triangular "edges" of a complete graph.

# Three-element subsets of six elements.
vertex.6 = t(combn(6,3))

x1 = dnorm(0:20, mean=10, sd=5)
x1 = (x1 - min(x1)) / max(x1)

# Transforms from 3-D to screen coordinates.
# x: the user coordinates, as a three-row matrix
#		(or a column vector), with rows:
#   x: in [-1,1]
#     FIXME: scale to Hasse diagram width?
#   y: in [0,20], as it's the number of triangles present
#   z: in [0,10]: this is arbitrary heights (for illustration;
#     they're presumably inaccurate)
p = function(x) {
	A = rbind(c(9, 0.3, 0), c(0, 0.2, 1))
	A %*% x
}

# Plots counts, sort of.
plot.counts = function() {
	b1 = rbind(x1, 0:20, 0)	
	p1 = t(p(b1))
	lines(p1, col="#00000040")
	b2 = rbind(-x1, 0:20, 0)	
	p2 = t(p(b2))
	lines(p2, col="#00000040")
	for(i in 1:21)
		lines(rbind(p1[i,], p2[i,]), col="#00000040")
}

# Plots a subset of triangles.
#   center: center of where to draw the triangles
#   i: row indices to include
#   r: radius of circle on which to put the points
#		label: label for the set
#   tri.colors: colors for the triangles (
#   Side effects: draws a subset of triangles
plot.tri.subset = function(center, i, r=0.8, label="",
		tri.colors=tri.colors) {
	n = 6
  circle.points = r * cbind(cos(2*pi*(1:n)/n), sin(2*pi*(1:n)/n))
#  print(circle.points)
	p = t(t(circle.points) + center)
  draw.triangles(p,
		vertex.6[i,,drop=FALSE],
    tri.colors$h[i], tri.colors$v[i])
	text(center[1]-r, center[2]+r, label)
}

# Wrapper for the above which:
#   - uses the projection transformation
#   - draws a line down to the "ground"
#   x, z: coordinates for the set
#   i: the set (which will determine the y coordinate)
#   label: the label for the set
plot.tri.subset.3d = function(x, z, i, label="", r=0.8) {
	# create the coordinates
	y = length(unique(i))
	p0 = as.vector(p(t(t(c(x, y, 0)))))
	p1 = as.vector(p(t(t(c(x, y, z)))))
	lines(rbind(p0, p1), col="#00000040", lwd=3)
	# for cheap perspective, make distant sets a bit smaller
	plot.tri.subset(p1, i, r=r*(1-y/70), label=label,
		tri.colors=tri.colors)
}

# Similar to plot.tri.subset.3d(), but plots triangles in
# some hue, with some of them greyed out
#   x0, x1, z: coordinates for the set
#   hue: the hue for the set
#   i: the larger set (as a vector of numbers)
#   j: the greyed-out subset (as a vector of numbers)
#     this will determine the y coordinate)
#   label: the label for the set
plot.tri.subset.3d.grey = function(x0, x1, z, hue, i, j, label="", r=0.8) {
	# make sure triangle indices in j are a subset of those
	# in i (if they aren't already)
	i = setdiff(i, j)
	# create the coordinates
	y = length(unique(j))
	p0 = as.vector(p(t(t(c(x0, y, 0)))))
	p1 = as.vector(p(t(t(c(x1, y, z)))))
	lines(rbind(p0, p1), col="#00000040", lwd=3)
	# first, the sets of whatever color (if any)
	if (length(i) > 0) {
		colors1 = data.frame(h=rep(hue,20), v=rep(0.7,20))
		# for cheap perspective, make distant sets a bit smaller
		plot.tri.subset(p1, i, r=r*(1-y/70), label=label,
			tri.colors=colors1)
	}
	# then, the grey sets (if any)
	if (length(j) > 0) {
		colors2 = data.frame(h=rep(0,20), v=rep(0,20))
		plot.tri.subset(p1, j, r=r*(1-y/70), label=label,
			tri.colors=colors2)
	}
}

# Hasse diagram, of some color-coded subsets.
pdf("Hasse.pdf", width=7, height=5)
par(mar=c(0,0,0,0))
world.bounds = cbind(c(-1,0,0), c(1,20,10))
# screen.bounds = range(p(world.bounds))

plot(0,0,
#	xlim=range(screen.bounds[1,]), ylim=range(screen.bounds[2,]),
#	xlim = screen.bounds, ylim = screen.bounds,
  xlim=c(-6,11), ylim=c(-0.2,12.2),
	type="n", xaxt="n", yaxt="n",xlab="", ylab="", bty="n")

plot.counts()

plot.tri.subset.3d(0, 0.6, c(2), "a)")
plot.tri.subset.3d(0, 2, c(1,2,5,11), "b)")
v = which(apply(vertex.6, 1, function(a) all(a!=6)))
plot.tri.subset.3d(0, 4.5, v, "c)")
plot.tri.subset.3d(0, 7, c(1:20), "d)")
plot.tri.subset.3d(-0.4, 2.7, c(1,2,5,11,3), "e)")
plot.tri.subset.3d(0.35, 2.9, c(1,2,3,4), "f)")
dev.off()

# Hasse diagram, showing omitting edges.
pdf("HasseWithOmissions.pdf", width=7, height=5)
par(mar=c(0,0,0,0))
world.bounds = cbind(c(-1,0,0), c(1,20,10))
plot(0,0,
  xlim=c(-6,11), ylim=c(-0.2,12.2),
	type="n", xaxt="n", yaxt="n",xlab="", ylab="", bty="n")
plot.counts()

# set of triangles, intermediate in size
s1 = c(1,2,3,4,5,6,11)
s2 = setdiff(s1, 4)

# at level 0
plot.tri.subset.3d.grey(0, -0.5, 3,
	0, c(1:20), c(), "a)")
plot.tri.subset.3d.grey(0, -0.2, 2.5,
	1/3, s1, c(), "")
plot.tri.subset.3d.grey(0, 0.19, 0.8,
	2/3, c(1:3), c(), "")

# at level 2
plot.tri.subset.3d.grey(0, -0.33, 5.5,
  0, c(1:20), c(1,2), "b)")
plot.tri.subset.3d.grey(0, -0.05, 2.5,
  1/3, s1, c(1,2), "")
plot.tri.subset.3d.grey(0, 0.33, 0.8,
  2/3, c(1:3), c(1,2), "")

# somewhere in the middle
plot.tri.subset.3d.grey(-0.1, 0, 5,
	0, c(1:20), s2, "c)")
plot.tri.subset.3d.grey(0.1, 0.1, 2.5,
	1/3, s1, s2, "")

# at level 19
plot.tri.subset.3d.grey(0, 0, 5,
	0, c(1:20), c(1:19), "d)")

# plot.tri.subset.3d.grey(0.1, 5, 2/3, c(1,2,3), c(1,2), "Z)")

# plot.tri.subset.3d.grey(0, 3, 0, c(1:20), c(), "Z)")
# plot.tri.subset.3d.grey(-0.1, 3, 0, c(1:20), c(1), "Z)")
# plot.tri.subset.3d.grey(-0.2, 3, 0, c(1:20),
#  	c(1,3,5,7,9,11,13,15,17), "Z)")

# plot.tri.subset.3d.grey(0, 0, 3, 0, c(1:20), c(20), "Z)")
dev.off()

