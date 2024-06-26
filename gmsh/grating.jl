"""
.geo script for grating

lc = 1e-1;

// Points
Point(1) = {0, 0, 0, lc};
Point(21) = {2*Pi, 0, 0, lc};
Point(22) = {2*Pi, Pi, 0, lc};
Point(23) = {0, Pi, 0, lc};

// use points to approximate a sine curve
h = 2*Pi/20.0;

For i In {2:20}
    Point(i) = {(i - 1)*h, Sin((i - 1)*h), 0, lc};
EndFor

// Lines
Line(1) = {23, 1};

For j In {2:23}
    Line(j) = {j - 1, j};
EndFor

Periodic Curve {22} = {1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
                 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 
                 23};

// Surface
Plane Surface(1) = {1};
"""

using Gmsh

gmsh.initialize()

gmsh.model.add("grating")
lc = 1e-1

gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(2*pi, 0, 0, lc, 21)
gmsh.model.geo.addPoint(2*pi, pi, 0, lc, 22)
gmsh.model.geo.addPoint(0, pi, 0, lc, 23)

# use Points to approximate a sine curve
h = 2*pi/20.0

for i in 2:20
    gmsh.model.geo.addPoint((i - 1)*h, sin((i - 1)*h), 0, lc, i)
end

gmsh.model.geo.addLine(23, 1, 1)

for i in 2:23
    gmsh.model.geo.addLine(i - 1, i, i)
end

gmsh.model.geo.addCurveLoop([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 
23], 1)

gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("grating.msh")
gmsh.finalize()
