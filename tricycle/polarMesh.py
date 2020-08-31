# Utilities for meshes using polar coordinates.

import math

import Part

# constructs a rectangular face
def quadFace(p):
  return [p[0], p[2], p[3], p[0], p[1], p[2]]

# converts to polar coordinates
def polar(x):
	th = x[0]
	r = x[1]
	z = x[2]
	return [ r * math.cos(th), r * math.sin(th), z ]

# constructs a rectangular face, using polar coordinates
def quadFacePolar(p):
	corners =[polar(x) for x in p]
	return quadFace(corners)

# returns one plane around the circle
# note float() cast, since FreeCAD is using Python 2
def facePolar(i):
	th1 = 2 * math.pi * (float(i) / numDivs)
	th2 = 2 * math.pi * (float(i+1) / numDivs)
	return quadFacePolar([
		[th1, 1, 1], [th1, 1, 0], [th2, 1, 0], [th2, 1, 1]])

# Constructs a surface by sweeping a shape (admittedly I
# could do this in FreeCAD, but this seems easier)
# Args:
#   r, z - coordinates of shape to sweep
#		(the shape is not automatically closed, but this
#		effect can be obtained by repeating the first
#		coordinates at the end)
#   numDivs - the number of divisions to include
# Returns: coordinates which can be passed to Mesh.mesh()
def sweep(r, z, numDivs):
	a = []
	for i in range(0, numDivs):
		for j in range(0, len(r)-1):
			th0 = 2 * math.pi * (float(i) / numDivs)
			th1 = 2 * math.pi * (float(i+1) / numDivs)
			f = quadFacePolar([
				[th0, r[j], z[j]],
				[th0, r[j+1], z[j+1]],
				[th1, r[j+1], z[j+1]],
				[th1, r[j], z[j]]
			])
			a.extend(f)
	return a

# Convert from a Mesh to a Shape (possibly not used)
def mesh2shape(m):
	shape = Part.Shape()
	shape.makeShapeFromMesh(m.Topology, 0.05)
	return Part.makeSolid(shape)

