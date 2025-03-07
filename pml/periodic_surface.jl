#----------------------------------------------------------
# Incident field
#----------------------------------------------------------

function uⁱ(y, k, θ)
    return exp(-im*k*cos(θ)*y)
end

#----------------------------------------------------------
# PML routines
#----------------------------------------------------------

function coord_transform(y)
    # PML parameters 
    σ = 20
    δ = 2
    # Outside the PML layer
    s = 1.0 + 0.0im
    # In the PML layer 
    if y > 1 
        s = 1.0 + im*σ*((y - 1)/δ)^2
    end

    return s
end

#----------------------------------------------------------
# FEM routines
#----------------------------------------------------------

function setup_grid(lc=0.5)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(π/10, 0.1, 0, lc)
    p3 = gmsh.model.geo.addPoint(π/5, 0, 0, lc)
    p4 = gmsh.model.geo.addPoint(3π/10, 0.1, 0, lc)
    p5 = gmsh.model.geo.addPoint(2π/5, 0, 0, lc)
    p6 = gmsh.model.geo.addPoint(π/2, 0.1, 0, lc)
    p7 = gmsh.model.geo.addPoint(3π/5, 0, 0, lc)
    p8 = gmsh.model.geo.addPoint(7π/10, 0.1, 0, lc)
    p9 = gmsh.model.geo.addPoint(4π/5, 0, 0, lc)
    p10 = gmsh.model.geo.addPoint(9π/10, 0.1, 0, lc)
    p11 = gmsh.model.geo.addPoint(π, 0, 0, lc)
    p12 = gmsh.model.geo.addPoint(π, 3, 0, lc)
    p13 = gmsh.model.geo.addPoint(0, 3, 0, lc)
  
    # Add the lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p5)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p9)
    l9 = gmsh.model.geo.addLine(p9, p10)
    l10 = gmsh.model.geo.addLine(p10, p11)
    l11 = gmsh.model.geo.addLine(p11, p12)
    l12 = gmsh.model.geo.addLine(p12, p13)
    l13 = gmsh.model.geo.addLine(p13, p1)

    # Add the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the Physical Groups
    gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [l11], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l12], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l13], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")

    # Set periodic mesh
    transform = [1, 0, 0, π, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [l11], [l13], transform)

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    # grid = mktempdir() do dir
    #     path = joinpath(dir, "mesh.msh")
    #     gmsh.write(path)
    #     togrid(path) 
    # end
    msh_path = joinpath(pwd(), "surface.msh")
    vtk_path = joinpath(pwd(), "surface.vtk")
    gmsh.write(msh_path)
    gmsh.write(vtk_path)
    grid = togrid(msh_path)

    # Finalize the Gmsh library
    gmsh.finalize()
    return grid
end

function setup_fevs(ip)
    qr = QuadratureRule{RefTriangle}(2)
    cv = CellValues(qr, ip)
    return cv
end

function setup_dofs(grid::Grid, ip)
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)
    return dh
end

# Be careful of the Dirichlet condition on "bottom" and the period π
function setup_bdcs(dh::DofHandler, k, θ)
    cst = ConstraintHandler(ComplexF64, dh)

    # Set periodic boundary condition on the "right" and "left"
    pfacets = collect_periodic_facets(dh.grid, "right", "left", x -> x + Vec{2}((π, 0.0)))
    pbc = PeriodicDirichlet(:u, pfacets)
    add!(cst, pbc)

    # Set Dirichlet boundary condition on the "bottom" and "top"
    dbc_top = Dirichlet(:u, getfacetset(dh.grid, "top"), x -> 0.0 + 0.0im)
    add!(cst, dbc_top)
    dbc_bottom = Dirichlet(:u, getfacetset(dh.grid, "bottom"), x -> -uⁱ(x[2], k, θ))
    # dbc_bottom = Dirichlet(:u, getfacetset(dh.grid, "bottom"), x -> 0.5)
    add!(cst, dbc_bottom)

    close!(cst)
    return cst
end

function allocate_matrices(dh::DofHandler, cst::ConstraintHandler)
    sp = init_sparsity_pattern(dh)
    add_cell_entries!(sp, dh)
    add_constraint_entries!(sp, cst)
    M = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return M
end

function assemble_global!(cv::CellValues, dh::DofHandler, A::SparseMatrixCSC, k, θ)
    α = k*sin(θ)
    # Allocate the local stiffness matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)

    # Create an assembler
    assembler = start_assemble(A)

    # Loop over all cells 
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        # Assemble local stiffness matrix
        assemble_local!(cv, cell, Ae, k, α)
        # Assemble Ae into A 
        assemble!(assembler, celldofs(cell), Ae)
    end

    return A
end

function assemble_local!(cv::CellValues, cache::CellCache, Ae::Matrix, k, α)
    n_basefuncs = getnbasefunctions(cv)
    # Reset local stiffness matrix to 0
    fill!(Ae, 0.0 + 0.0im)
    coords = getcoordinates(cache)

    # Loop over quadrature points 
    for qp in 1:getnquadpoints(cv)
        dx = getdetJdV(cv, qp)
        coord_qp = spatial_coordinate(cv, qp, coords)
        s = coord_transform(coord_qp[2])
        # Loop over test shape functions
        for i in 1:n_basefuncs
            v = shape_value(cv, qp, i)
            ∇v = shape_gradient(cv, qp, i)
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                u = shape_value(cv, qp, j)
                ∇u = shape_gradient(cv, qp, j)
                # Assemble local stiffness matrix
                Ae[i, j] += (s*∇u[1]*∇v[1] + ∇u[2]*∇v[2]/s - 2im*α*s*∇u[1]*v - (k^2 - α^2)*s*u*v) * dx
            end
        end
    end

    return Ae
end


