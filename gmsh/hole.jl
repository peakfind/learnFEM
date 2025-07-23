using Gmsh

period = 2π
ratio = 0.2

# Initialize gmsh 
gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 2)

# Get the radius
r = period * ratio

# Add the points
# For the square
p1 = gmsh.model.geo.addPoint(period/2, -period/2, 0, lc)
p2 = gmsh.model.geo.addPoint(period/2, period/2, 0, lc)
p3 = gmsh.model.geo.addPoint(-period/2, period/2, 0, lc)
p4 = gmsh.model.geo.addPoint(-period/2, -period/2, 0, lc)

# For the circle
p5 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p6 = gmsh.model.geo.addPoint(r, 0, 0, lc)
p7 = gmsh.model.geo.addPoint(0, r, 0, lc)
p8 = gmsh.model.geo.addPoint(-r, 0, 0, lc)
p9 = gmsh.model.geo.addPoint(0, -r, 0, lc)

# Add the lines
# For the square
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# For the circle
c1 = gmsh.model.geo.addCircleArc(p6, p5, p7)
c2 = gmsh.model.geo.addCircleArc(p7, p5, p8)
c3 = gmsh.model.geo.addCircleArc(p8, p5, p9)
c4 = gmsh.model.geo.addCircleArc(p9, p5, p6)

# Create the loops and the domain
square = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
circle = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])
domain = gmsh.model.geo.addPlaneSurface([square, circle])

# Synchronize the model
gmsh.model.geo.synchronize()

# Generate a 2D mesh
gmsh.model.mesh.generate(2)

gmsh.write("column.vtk")
gmsh.finalize()