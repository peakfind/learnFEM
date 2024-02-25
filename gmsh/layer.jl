using gmsh_jll
include(gmsh_jll.gmsh_api)

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

#TODO: how to set periodic boundary condtions?