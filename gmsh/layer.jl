"""
.geo script for layer

                           Line(4)
     P4(0,pi) ------------------------------- P3(2π,π)
              |                             |
              |                             |
   Line(1)    |                             |  Line(3)
              |                             |
              |                             |
              |                             |
      P1(0,0) ------------------------------- P2(2π,0)
                           Line(2)
// Elementary entity: point
// A Point is uniquely identified by:
//     a tag: a strictly positive integer.
//     three coordinates: X, Y and Z.
//     lc: the target mesh size close to the point.

lc = 1e-1;
Point(1) = {0, 0, 0, lc};
Point(2) = {2*Pi, 0, 0, lc};
Point(3) = {2*Pi, Pi, 0, lc};
Point(4) = {0, Pi, 0, lc};

// The second elementary entity: curve
// A straight line is identified by:
//     a tag: curve tags are separate from point tags.
//     a list of two point tags: start point -> end point.

Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

// The third elementary entity: surface
// In order to define a simple rectangular surface from the four 
// curves defined above, a curve loop has first to be defined.
//
// A curve loop is identified by:
// a tag: unique amongst curve loops.
// an ordered list of connected curves: a sign being associated with 
// each curve (depending on the orientation of the curve to form a loop).

Curve Loop(1) = {1, 2, 3, 4};

// We can then define the surface as a list of curve loops

Plane Surface(1) = {1};

// Now Gmsh knows everything to display the geometry and to mesh it.
// An optional step is that we can group elementary geometrical entities 
// into more meaningful groups:
//     mathematical: domain, boundary.
//     functional: left wing, fuselage.
//     material: steel, carbon.
//
// Such groups in Gmsh are called "Physical Groups"
"""

using Gmsh

gmsh.initialize()

gmsh.model.add("layer")
lc = 0.5

gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(2*pi, 0, 0, lc, 2)
gmsh.model.geo.addPoint(2*pi, pi, 0, lc, 3)
gmsh.model.geo.addPoint(0, pi, 0, lc, 4)

gmsh.model.geo.addLine(4, 1, 1)
gmsh.model.geo.addLine(1, 2, 2)
gmsh.model.geo.addLine(2, 3, 3)
gmsh.model.geo.addLine(3, 4, 4)

gmsh.model.geo.addCurveLoop([1, 2, 3, 4],1)

gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("layer.msh")
gmsh.finalize()