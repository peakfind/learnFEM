# ---------------------------------------------------------
# pml.jl
# 
# Use PML to truncate the open waveguide problem to a closed 
# waveguide.
#
# Reference:
# A. Kirsch and R. Zhang, Computation of the exceptional values 
# for an open waveguide, preprint. Example 2 in Section 5.4.1
# ---------------------------------------------------------

using Gmsh
using Ferrite, FerriteGmsh
using SparseArrays
using Arpack
using Plots

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


"""
    setup_grid(;d=2π, h=1.5, δ=0.5, lc=0.5)

Generate a mesh by julia interface for Gmsh and read the `.msh` file by `FerriteGmsh`

# Arguments

- `d`: the period of the periodic layer
- `h₀`: the height of the periodic layer
- `δ`: the height of the PML layer
- `lc`: target mesh size close to the point

# Points

```
 4  -------------  3
 |                 |
 1  -------------  2
```

# Lines

```
 .  ----  3  ----  .
 4                 2 
 .  ----  1  ----  .
```

# Misc 

- We can generate a `*.vtk` file by replace the `do end` block with `gmsh.write(*.vtk)`
"""
function setup_grid(;d=2π, h=1.5, δ=0.5, lc=0.5)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Compute the height of the truncated domain
    b = h + δ

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(d, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(d, b, 0, lc) 
    p4 = gmsh.model.geo.addPoint(0, b, 0, lc)

    # Add the lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Add the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the Physical Groups
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [l2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")

    # Set periodic mesh
    transform = [1, 0, 0, d, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], transform)

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    grid = mktempdir() do dir
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path) 
    end

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



"""
    get_dispersion(cv::CellValues, dof::DofHandler, cst::ConstraintHandler, A::SparseMatrixCSC, B::SparseMatrixCSC, qms, wh::Symbol; neigs=3)

# Arguments

- `qms`: discreted quasi momentum
- `neigs`: the number of the expected eigenvalues

"""
function get_dispersion(cv::CellValues, dof::DofHandler, cst::ConstraintHandler, A::SparseMatrixCSC, B::SparseMatrixCSC, qms; wh::Symbol, neigs=3)
    # preallocate for the dispersion curves
    m = length(qms)
    dispersion = zeros(ComplexF64, m, neigs)

    # Assemble the matrix B
    B = assemble_b(cv, dof, B)
    apply!(B, cst)

    # For fixed α
    for (i, qm) in enumerate(qms)
        # Assemble the matrix A(α)
        A = assemble_a(cv, dof, A, qm)
        apply!(A, cst)

        # solve the general eigenvalue problem by Arpack
        λ, _ = eigs(A, B; nev=neigs, which = wh, maxiter = 3000) 

        dispersion[i, :] = λ
        fill!(A, 0.0 + 0.0im)
    end
    
    return dispersion
end

function assemble_a(cv::CellValues, dof::DofHandler, A::SparseMatrixCSC, qm)
    ndofs_c = ndofs_per_cell(dof)
    assembler = start_assemble(A)
    
    # Preallocate the local stiffness matrix
    Ae = zeros(ComplexF64, ndofs_c, ndofs_c)

    for cell in CellIterator(dof)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        fill!(Ae, 0)
        coords = getcoordinates(cell)

        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            s = coord_transform(coords_qp[2])
            for i in 1:ndofs_c
                v = shape_value(cv, qp, i)
                grad_v = shape_gradient(cv, qp, i)
                for j in 1:ndofs_c
                    u = shape_value(cv, qp, j) 
                    grad_u = shape_gradient(cv, qp, j)
                    Ae[i, j] += ((s*grad_u[1]*grad_v[1] + grad_u[2]*grad_v[2]/s) - (2*im*qm*s*grad_u[1]*v) + (qm^2*s*u*v)) * dx
                end
            end
        end
        
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A
end

function assemble_b(cv::CellValues, dof::DofHandler, B::SparseMatrixCSC)
    ndofs_c = ndofs_per_cell(dof)
    assembler = start_assemble(B)
    
    # Preallocate the local stiffness matrix
    Be = zeros(ComplexF64, ndofs_c, ndofs_c)

    for cell in CellIterator(dof)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        fill!(Be, 0)
        coords = getcoordinates(cell)
        
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            s = coord_transform(coords_qp[2])
            n = refractive_index(coords_qp[2])
            for i in 1:ndofs_c
                v = shape_value(cv, qp, i)
                for j in 1:ndofs_c
                    u = shape_value(cv, qp, j)
                    Be[i, j] += (n * s * u * v) * dx
                end
            end
        end
        assemble!(assembler, celldofs(cell), Be)
    end

    return B
end

function main()
    d = 2π
    grid = setup_grid(;d=2π, h=1.5, δ=0.5, lc=0.5)

    # Define the interpolation: linear lagrange
    ip = Lagrange{RefTriangle, 1}()

    cv = setup_fevs(ip)
    dh = setup_dofs(grid, ip)
    ch = setup_bdcs(dh, d)

    # Allocate the matrices
    A = allocate_matrices(dh, ch)
    B = allocate_matrices(dh, ch)

    # Brillouin zone α
    α = collect(range(-pi/d, pi/d, 40))

    # Solve the eigenvalue problem
    k = get_dispersion(cv, dh, ch, A, B, α; wh=:SM, neigs=3)
    
    @show k
end

main()