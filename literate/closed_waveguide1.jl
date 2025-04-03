# # Dispersion diagram of closed waveguides I

# ## Introduction 

# ## Implementation

using Gmsh, Ferrite, FerriteGmsh, SparseArrays, Arpack, Plots

# ### Mesh generation with `Gmsh.jl` 
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

# ## Plain program

# Here follows a version of the program without any comments.
# ```julia
# using Gmsh, Ferrite, FerriteGmsh, SparseArrays, Arpack, Plots
# 
# function setup_grid(lc = 0.05)
#     gmsh.initialize()
#     gmsh.option.setNumber("General.Verbosity", 2)
#     
#     p1 = gmsh.model.geo.addPoint(-1/2, 0, 0, lc)
#     p2 = gmsh.model.geo.addPoint(1/2, 0, 0, lc)
#     p3 = gmsh.model.geo.addPoint(1/2, 1, 0, lc)
#     p4 = gmsh.model.geo.addPoint(-1/2, 1, 0, lc)
#     
#     l1 = gmsh.model.geo.addLine(p1, p2)
#     l2 = gmsh.model.geo.addLine(p2, p3)
#     l3 = gmsh.model.geo.addLine(p3, p4)
#     l4 = gmsh.model.geo.addLine(p4, p1)
#     
#     loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
#     surf = gmsh.model.geo.addPlaneSurface([loop])
#     
#     gmsh.model.geo.synchronize()
#     
#     gmsh.model.addPhysicalGroup(1, [l1], -1, "Gamma")
#     gmsh.model.addPhysicalGroup(1, [l2], -1, "Gamma1")
#     gmsh.model.addPhysicalGroup(1, [l3], -1, "Gammab")
#     gmsh.model.addPhysicalGroup(1, [l4], -1, "Gamma2")
#     gmsh.model.addPhysicalGroup(2, [surf], -1, "Omegab")
#     
#     gmsh.model.mesh.setPeriodic(1, [l2], [l4], [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]) 
#     
#     gmsh.model.mesh.generate(2)
#     
#     path = joinpath(pwd(), "cell.msh")
#     gmsh.write(path)
#     grid = togrid(path)
#     
#     gmsh.finalize()
#     
#     return grid
# end
# 
# function setup_vals(interpolation)
#     qr = QuadratureRule{RefTriangle}(2)
#     cellvalues = CellValues(qr, interpolation)
#     return cellvalues
# end
# 
# function setup_dofs(grid::Grid, interpolation)
#     dh = DofHandler(grid)
#     add!(dh, :u, interpolation)
#     close!(dh)
#     return dh
# end
# 
# function setup_bcs(dofhandler::DofHandler)
#     ch = ConstraintHandler(dofhandler)
#     periodic_faces = collect_periodic_facets(dofhandler.grid, "Gamma1", "Gamma2", x -> x + Vec{2}((1.0, 0.0)))
#     pbc = PeriodicDirichlet(:u, periodic_faces)
#     add!(ch, pbc)
#     close!(ch)
#     return ch
# end
# 
# function allocate_matries(dofhandler::DofHandler, csthandler::ConstraintHandler)
#     sp = init_sparsity_pattern(dofhandler)
#     add_cell_entries!(sp, dofhandler)
#     add_constraint_entries!(sp, csthandler)
#     K = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
#     return K
# end
# 
# function assemble_a(cellvalues::CellValues, A::SparseMatrixCSC, dh::DofHandler, α)
#     ndofs_c = ndofs_per_cell(dh)
#     assembler = start_assemble(A)
# 
#     Ae = zeros(ComplexF64, ndofs_c, ndofs_c)
# 
#     for cell in CellIterator(dh)
#         reinit!(cellvalues, cell)
#         fill!(Ae, 0)
#         
#         for qp in 1:getnquadpoints(cellvalues)
#             dx = getdetJdV(cellvalues, qp)
#             for i in 1:ndofs_c
#                 v = shape_value(cellvalues, qp, i)
#                 grad_v = shape_gradient(cellvalues, qp, i)
#                 for j in 1:ndofs_c
#                     u = shape_value(cellvalues, qp, j) 
#                     grad_u = shape_gradient(cellvalues, qp, j)
#                     Ae[i, j] += (dot(grad_u, grad_v) - (2*im*α*grad_u[1]*v) + (α^2*u*v)) * dx
#                 end
#             end
#         end
#         
#         assemble!(assembler, celldofs(cell), Ae)
#     end
# 
#     return A
# end
# 
# function assemble_b(cellvalues::CellValues, B::SparseMatrixCSC, dh::DofHandler)
#     ndofs_c = ndofs_per_cell(dh)
#     assembler = start_assemble(B)
# 
#     Be = zeros(ComplexF64, ndofs_c, ndofs_c)
# 
#     for cell in CellIterator(dh)
#         reinit!(cellvalues, cell)
#         fill!(Be, 0)
#         
#         for qp in 1:getnquadpoints(cellvalues)
#             dx = getdetJdV(cellvalues, qp)
#             for i in 1:ndofs_c
#                 v = shape_value(cellvalues, qp, i)
#                 for j in 1:ndofs_c
#                     u = shape_value(cellvalues, qp, j)
#                     Be[i, j] += (u*v) * dx
#                 end
#             end
#         end
#         assemble!(assembler, celldofs(cell), Be)
#     end
# 
#     return B
# end
# 
# function cal_eigs(cellvalues::CellValues, dh::DofHandler, csthandler::ConstraintHandler, 
#               A::SparseMatrixCSC, B::SparseMatrixCSC, f, brillouin_zone, nevs)
#     m = length(brillouin_zone)
#     disperion_diagram = zeros(m, nevs)
#     B = assemble_b(cellvalues, B, dh)
#     apply!(B, f, csthandler)
# 
#     for (ind, α) in enumerate(brillouin_zone)
#         A = assemble_a(cellvalues, A, dh, α) 
#         apply!(A, f, csthandler)
#         
#         vals, _ = eigs(A, B, nev = nevs, which = :SR)
#         vals = real(vals)
#         disperion_diagram[ind,:] = vals
#         fill!(A, 0)
#     end
# 
#     return disperion_diagram
# end
# 
# function main()
#     grid = setup_grid(0.04)
#     ip = Lagrange{RefTriangle, 1}()
#     cv = setup_vals(ip)
#     dh = setup_dofs(grid, ip)
#     ch = setup_bcs(dh)
#     
#     A = allocate_matries(dh, ch)
#     B = allocate_matries(dh, ch)
#     f = zeros(ComplexF64, ndofs(dh))
#     
#     brillouin = collect(range(-pi, pi, 100))
#     disp_diag = cal_eigs(cv, dh, ch, A, B, f, brillouin, 8)
#     plot(brillouin, disp_diag, 
#          legend = false,
#          xlims = (-pi, pi),
#          ylims = (0,80),
#          linewidth = 2,
#          title="Dispersion diagram")
# end
# 
# main()
# ```