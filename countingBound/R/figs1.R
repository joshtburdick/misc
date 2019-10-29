# Draws some figures.

# default triangle colors
set.seed(1)
tri.colors = data.frame(h=runif(20), v=runif(20, min=0.6,1))

# Draws some triangles.
#   x: matrix with each row as coordinates
#   corners: indices to use (with each row being three row
#     indices into x)
#   h, v: hue and value of colors to use
draw.triangles = function(x, corners, h, v) {
  border = hsv(h, 1, v, 0.6)
  fill = hsv(h, 0.5, v, 0.2)
  for(i in 1:nrow(corners)) {
    polygon(x[ corners[i,] , ], border = border[i], col = fill[i], lwd=2)  
  }
}

# Plots triangular "edges" of a complete graph.
plot.k = function(n) {
  lim = c(-1.2,1.2)
  plot(0,0, xlim=lim, ylim=lim, type="n",
    xaxt="n", yaxt="n")
  circle.points = cbind(cos(2*pi*(1:n)/n), sin(2*pi*(1:n)/n))
  draw.triangles(circle.points, t(combn(n,3)),
    tri.colors$h, tri.colors$v)
}

# Plots some "less overlapping" triangles.
plot.less.overlapping = function() {
  p = cbind(x=c(-0.3,0.3, rnorm(20)),
    y=c(-3, -3, rnorm(20)))

  corners = cbind(1, 2, 3:22)
  plot(0,0, xlim=c(-2.5, 2.5), ylim=c(-3,2),
    type="n", xaxt="n", yaxt="n")
  
  draw.triangles(p, corners, tri.colors$h, tri.colors$v)
}

# Plots some partially-overlapping cliques.
plot.red.blue = function(offset) {
  n = 8
  circle.points = cbind(cos(2*pi*(1:n)/n), sin(2*pi*(1:n)/n))
  lim = c(-1.2,1.2)
  plot(0,0, xlim=lim, ylim=lim, type="n",
    xaxt="n", yaxt="n")
  draw.triangles(circle.points,
    t(combn(4,3)),
    rep(0,4), rep(0.5,5))
  draw.triangles(circle.points,
    t(combn(4,3)) + offset,
    rep(2/3,4), rep(0.8,5))
}

pdf("tri1.pdf", width=7.5, height=1.5)
layout(cbind(1,2), width=c(1,4))
par(mar=c(0,0,0,0), bty="n")
plot.k(6)
plot.less.overlapping()
dev.off()

# partially-overlapping cliques
pdf("overlapping.pdf", width=6, height=1.5)
par(mfrow=c(1,4), mar=c(0,0,0,0), bty="n")
for(i in 1:4)
  plot.red.blue(i)
dev.off()

# example of distinction between maximal and optimal coverings
pdf("maximal.pdf", width=5, height=2)
par(mfrow=c(1,2), mar=c(0,0,0,0)+0.2, bty="n")
x = cbind(c(-2:2,0), c(1,0,1,0,1,-1)*sqrt(3)) / 2
x = rbind(x, x+0.04)  # points, slightly offset
# first, plot all three triangles
plot(0,0, xlim=range(x[,1]), ylim=range(x[,2]), type="n", xaxt="n", yaxt="n")
draw.triangles(x, rbind(c(1,2,3), c(2,4,6), c(3,4,5), c(8,9,10)),
  c(0,1,2,4)/6, c(0.5,0.5,0.5,0.8))
# then, omit the triangle in the center
plot(0,0, xlim=range(x[,1]), ylim=range(x[,2]), type="n", xaxt="n", yaxt="n")
draw.triangles(x, rbind(c(1,2,3), c(2,4,6), c(3,4,5)),
  c(0,1,2,4)/6, c(0.5,0.5,0.5,0.4))
dev.off()

# example of the aforementioned bound
pdf("mesh2D.pdf", width=5, height=5)
par(mar=c(1,1,1,1))
plot(0,0, xlim=c(0,10), ylim=c(0,10), type="n",	xaxt="n", yaxt="n")
for (i in 1:10)
	lines(c(i,0), c(0,11-i))
abline(0,1, col="red")
# FIXME shade the feasible regions?
dev.off()

# various Hasse diagrams
source("Hasse.R")

