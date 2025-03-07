# The coordinate transformation
function coord_transform(x; h=1.5, δ=0.5, χ=1.0+2.0im, m=2, ρ=80)
    if x <= h
        return 1.0 + 0.0im
    else 
        return 1.0 + ρ*χ*(((x - h)/δ)^m)
    end
end

# The refractive index
function refractive_index(x; h₀=1.0)
    if x < h₀
        return 3.9
    else
        return 1.0
    end
end

# Setup the grid
function setup_grid(;d=2π, h=1.5, δ=0.5, lc=0.5)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    b = h + δ
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(d, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(d, b, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, b, 0, lc)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [l2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")
    transform = [1, 0, 0, d, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], transform)
    gmsh.model.mesh.generate(2)
    grid = mktempdir() do dir
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path) 
    end
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

function setup_bdcs(dof::DofHandler, d)
    cst = ConstraintHandler(dof)

    # Set periodic boundary condition on the "right" and "left"
    pfacets = collect_periodic_facets(dof.grid, "right", "left", x -> x + Vec{2}((d, 0.0)))
    pbc = PeriodicDirichlet(:u, pfacets)
    add!(cst, pbc)

    # Set Dirichlet boundary condition on the "bottom" and "top"
    dfacets = union(getfacetset(dof.grid, "bottom"), getfacetset(dof.grid, "top"))
    dbc = Dirichlet(:u, dfacets, x -> 0)
    add!(cst, dbc)

    close!(cst)
    return cst
end

function allocate_matrices(dof::DofHandler, cst::ConstraintHandler)
    sp = init_sparsity_pattern(dof)
    add_cell_entries!(sp, dof)
    add_constraint_entries!(sp, cst)
    M = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return M
end

function assemble_A2(cv::CellValues, dof::DofHandler, A₂::SparseMatrixCSC)
    ndofs_c = ndofs_per_cell(dof)
    assembler = start_assemble(A₂)
    # Preallocate the local stiffness matrix
    Loc = zeros(ComplexF64, ndofs_c, ndofs_c)
    for cell in CellIterator(dof)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        fill!(Loc, 0.0 + 0.0im)
        coords = getcoordinates(cell)
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            s = coord_transform(coords_qp[2])
            for i in 1:ndofs_c
                v = shape_value(cv, qp, i)
                for j in 1:ndofs_c
                    u = shape_value(cv, qp, j)
                    Loc[i, j] += (s * u * v) * dx
                end
            end
        end
        assemble!(assembler, celldofs(cell), Loc)
    end
    return A₂
end

function assemble_A1(cv::CellValues, dof::DofHandler, A₁::SparseMatrixCSC)
    ndofs_c = ndofs_per_cell(dof)
    assembler = start_assemble(A₁)
    
    # Preallocate the local stiffness matrix
    Loc = zeros(ComplexF64, ndofs_c, ndofs_c)

    for cell in CellIterator(dof)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        fill!(Loc, 0)
        coords = getcoordinates(cell)

        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            s = coord_transform(coords_qp[2])
            for i in 1:ndofs_c
                v = shape_value(cv, qp, i)
                for j in 1:ndofs_c 
                    grad_u = shape_gradient(cv, qp, j)
                    Loc[i, j] += (-2im * s * grad_u[1] * v) * dx
                end
            end
        end
        
        assemble!(assembler, celldofs(cell), Loc)
    end
    
    return A₁
end

function assemble_A0(cv::CellValues, dof::DofHandler, A₀::SparseMatrixCSC)
    ndofs_c = ndofs_per_cell(dof)
    assembler = start_assemble(A₀)
    
    # Preallocate the local stiffness matrix
    Loc = zeros(ComplexF64, ndofs_c, ndofs_c)

    for cell in CellIterator(dof)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        fill!(Loc, 0)
        coords = getcoordinates(cell)

        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            s = coord_transform(coords_qp[2])
            n = refractive_index(coords_qp[2])
            for i in 1:ndofs_c
                v = shape_value(cv, qp, i)
                grad_v = shape_gradient(cv, qp, i)
                for j in 1:ndofs_c
                    u = shape_value(cv, qp, j) 
                    grad_u = shape_gradient(cv, qp, j)
                    Loc[i, j] += ((s * grad_u[1] * grad_v[1]) + (grad_u[2] * grad_v[2])/s - (4.1^2) * n * u * v) * dx
                end
            end
        end
        
        assemble!(assembler, celldofs(cell), Loc)
    end
    
    return A₀
end

function is_real(λ; tol=1e-12)
    if abs(imag(λ)) < tol
        return true
    else
        return false
    end
end