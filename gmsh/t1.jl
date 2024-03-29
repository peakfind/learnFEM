using gmsh_jll
include(gmsh_jll.gmsh_api)

gmsh.initialize()

gmsh.model.add("t1")
lc = 1e-2

gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(.1, 0, 0, lc, 2)
gmsh.model.geo.addPoint(.1, .3, 0, lc, 3)
gmsh.model.geo.addPoint(0, .3, 0, lc, 4)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(3, 2, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.addPhysicalGroup(0, [1, 2], 1)
gmsh.model.addPhysicalGroup(1, [1, 2], 2)
gmsh.model.addPhysicalGroup(2, [1], 6)
gmsh.model.setPhysicalName(2, 6, "My surface")

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("t1.msh")
gmsh.finalize()