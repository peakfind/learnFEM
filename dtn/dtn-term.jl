#=
# dtn-term.jl 
#
# Assemble the dtn terms and add it to the global stiffness matrix
=#

using Ferrite, FerriteGmsh, Gmsh

function setup_grid(lc=0.5)
    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(2, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(2, 2, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, 2, 0, lc)

    # Add the lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Create the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the physical domains (for boundary conditions)
    gmsh.model.addPhysicalGroup(1, [l1], -1, "Gamma") 
    gmsh.model.addPhysicalGroup(1, [l2], -1, "Gamma1")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "Gammab")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "Gamma2")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Omegab")

    # Set Periodic boundary condition
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], [1, 0, 0, 2, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    # Save the mesh and read the .msh file in Ferrite
    path = joinpath(pwd(), "mesh.msh")
    path_ = joinpath(pwd(), "mesh.vtk")
    gmsh.write(path)
    gmsh.write(path_)
    grid = togrid(path)

    # Finalize the Gmsh library
    gmsh.finalize()

    return grid
end

grid = setup_grid()