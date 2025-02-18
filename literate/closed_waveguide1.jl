# # Dispersion diagram of closed waveguides I
#
#-
#md # !!! tip
#md #     The code in this note is written with `Ferrite.jl` v1.0.
#-
#
# ## Introduction and problem formulation
# 
# 
# ## Commented program
#

using Gmsh, Ferrite, FerriteGmsh, SparseArrays, Arpack, Plots

# ### Mesh generation with `Gmsh.jl`
# 
# 

function setup_grid(lc = 0.05)
    ## Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    
    ## Add the points
    p1 = gmsh.model.geo.addPoint(-1/2, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(1/2, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(1/2, 1, 0, lc)
    p4 = gmsh.model.geo.addPoint(-1/2, 1, 0, lc)

    ## Add the lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    ## Create the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    ## Synchronize the model
    gmsh.model.geo.synchronize()

    ## Create the physical domains
    gmsh.model.addPhysicalGroup(1, [l1], -1, "Gamma")
    gmsh.model.addPhysicalGroup(1, [l2], -1, "Gamma1")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "Gammab")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "Gamma2")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Omegab")

    ## Set Periodic boundary condition
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]) 

    ## Generate a 2D mesh
    gmsh.model.mesh.generate(2)
    
    ## Save the mesh and read the .msh file in Ferrite
    path = joinpath(pwd(), "cell.msh")
    gmsh.write(path)
    grid = togrid(path)
    
    ## Finalize the Gmsh library
    gmsh.finalize()
    
    return grid
end

# ### FEvalues

function setup_vals(interpolation)
    qr = QuadratureRule{RefTriangle}(2)
    cellvalues = CellValues(qr, interpolation)
    return cellvalues
end

# ### Degrees of freedom
function setup_dofs(grid::Grid, interpolation)
    dh = DofHandler(grid)
    add!(dh, :u, interpolation)
    close!(dh)
    return dh
end

# ### Periodic boundary conditions
function setup_bcs(dofhandler::DofHandler)
    ch = ConstraintHandler(dofhandler)
    ## Periodic boundary condition
    periodic_faces = collect_periodic_facets(dofhandler.grid, "Gamma1", "Gamma2", x -> x + Vec{2}((1.0, 0.0)))
    pbc = PeriodicDirichlet(:u, periodic_faces)
    add!(ch, pbc)
    close!(ch)
    return ch
end

# ### Allocate Sparse matrix
function allocate_matries(dofhandler::DofHandler, csthandler::ConstraintHandler)
    sp = init_sparsity_pattern(dofhandler)
    add_cell_entries!(sp, dofhandler)
    add_constraint_entries!(sp, csthandler)
    K = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return K
end

# ### Assemble matrices