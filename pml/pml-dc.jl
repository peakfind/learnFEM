#=
# pml-dc.jl
#
# Use PML to truncate the unbounded domain.
# Compute the eigenvalues of smallest magnitude.
# 
# Author: Jiayi Zhang
# Date:   21/01/2025
=#

using Ferrite, FerriteGmsh
using SparseArrays, StaticArrays, Arpack
using Gmsh
using Plots

function coord_transform(x)
    # PML parameters
    b = 2.0
    δ = 0.5
    ρ = 2.0
    χ = exp(im*pi/4)
    m = 1.0

    if x < b
        return 1.0 + 0.0im
    else
        return 1.0 + ρ * χ * ((x - b)/δ)^m
    end
end

"""
    setup_grid(d, n, b, lc)

Generate the mesh by julia interface for Gmsh and read the `.msh` file by `FerriteGmsh`

# Arguments

- `d`: the period of the periodic grating
- `n`: the number of points on the profile of the grating
- `b`: the height of th truncated domain
- `lc`: target mesh size close to the point

# Points

```
n+2 ------------- n+1
 |                 |
 1  -------------  n
```

# Lines

```
 .  ---- n+1 ----  .
n+2                n 
 .  -- 1 ~ n-1 --  .
```

# Examples

```julia
d=2*pi, n=10, b=2
setup_grid(2*pi, 10, 2)
```

# Misc 

- We can generate a `*.vtk` file by replace the `do end` block with `gmsh.write(*.vtk)`
"""
function setup_grid(d, n, b, lc=0.5)
    h = d/(n-1)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    # Use a polyline to approximate the profile
    points = zeros(Int32, n+2)

    points[1] = gmsh.model.geo.addPoint(0, 0, 0, lc)
    points[n] = gmsh.model.geo.addPoint(d, 0, 0, lc)

    for i in 2:n-1
        points[i] = gmsh.model.geo.addPoint((i-1)*h, sin((i-1)*h), 0, lc) 
    end

    points[n+1] = gmsh.model.geo.addPoint(d, b, 0, lc) 
    points[n+2] = gmsh.model.geo.addPoint(0, b, 0, lc)

    # Add the lines
    lines = zeros(Int32, n+2)

    for i in 1:n+1
        lines[i] = gmsh.model.geo.addLine(points[i], points[i+1])
    end

    lines[n+2] = gmsh.model.geo.addLine(points[n+2], points[1])

    # Add the loop and the surface
    loop = gmsh.model.geo.addCurveLoop(lines)
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the Physical Groups
    gmsh.model.addPhysicalGroup(1, [lines[n+1]], -1, "Top")
    gmsh.model.addPhysicalGroup(1, [lines[n+2]], -1, "Left")
    gmsh.model.addPhysicalGroup(1, [lines[n]], -1, "Right")
    gmsh.model.addPhysicalGroup(1, lines[1:(n-1)], -1, "bottom")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Omega")

    # Set periodic mesh
    transform = @SVector [1, 0, 0, d, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [lines[n]], [lines[n+2]], transform)

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

function setup_fevs(interpolation)
    qr = QuadratureRule{RefTriangle}(2)
    cellvalues = CellValues(qr, interpolation)
    return cellvalues
end

function setup_dofs(grid::Grid, interpolation)
    dh = DofHandler(grid)
    add!(dh, :u, interpolation)
    close!(dh)
    return dh
end

function setup_bdcs(dofhandler::DofHandler, d)
    csthandler = ConstraintHandler(dofhandler)

    # Set periodic boundary condition on the "Right" and "Left"
    periodic_faces = collect_periodic_facets(dofhandler.grid, "Right", "Left", x -> x + Vec{2}((d, 0.0)))
    pbc = PeriodicDirichlet(:u, periodic_faces)
    add!(csthandler, pbc)

    # Set Dirichlet boundary condition on the "Top"
    dbc = Dirichlet(:u, getfacetset(dofhandler.grid, "Top"), x -> 0)
    add!(csthandler, dbc)

    close!(csthandler)
    return csthandler
end

function allocate_matrices(dofhandler::DofHandler, csthandler::ConstraintHandler)
    sp = init_sparsity_pattern(dofhandler)
    add_cell_entries!(sp, dofhandler)
    add_constraint_entries!(sp, csthandler)
    M = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return M
end

"""
    get_dispersion(cellvalues, dofhandler, csthandler, A, B, qm, neigs)

Formulate the general eigenvalue problem (GEP) 
```
    A(α)x = k^2 Bx
``` 
by FEM and solve it by `Arpack`

# Arguments

- `A`: the matrix `A(α)` in the GEP
- `B`: the matrix `B` in the GEP
- `qms`: the quasimomentums
- `neigs`: the number of the eigenvalues we want to compute (needed by `Arpack`)

# Steps

1. Fix `α` (entries in `qms`)
2. Assemble matrices `A` (resp. `B`) by the function `assemble_a()` (resp. `assemble_b`)
3. solve the eigenvalue problem
"""
function get_dispersion(cellvalues::CellValues, dofhandler::DofHandler, csthandler::ConstraintHandler, A::SparseMatrixCSC, B::SparseMatrixCSC, qms, neigs)
    m = length(qms)
    dispersion = zeros(ComplexF64, m, neigs)

    # Assemble the matrix B
    B = assemble_b(cellvalues, B, dofhandler)
    apply!(B, csthandler)

    # For fixed α
    for (i, qm) in enumerate(qms)
        # Assemble the matrix A(α)
        A = assemble_a(cellvalues, A, dofhandler, qm)
        apply!(A, csthandler)

        # solve the general eigenvalue problem
        λ, _ = eigs(A, B, nev = neigs, which = :SM)
        dispersion[i, :] = λ
        fill!(A, 0.0)
    end
    
    return dispersion
end

function assemble_a(cellvalues::CellValues, A::SparseMatrixCSC, dofhandler::DofHandler, qm)
    ndofs_c = ndofs_per_cell(dofhandler)
    assembler = start_assemble(A)
    coef = @MMatrix zeros(ComplexF64, 2, 2)
    
    # Preallocate the local stiffness matrix
    Ae = zeros(ComplexF64, ndofs_c, ndofs_c)

    for cell in CellIterator(dofhandler)
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        fill!(Ae, 0)
        coords = getcoordinates(cell)

        for qp in 1:getnquadpoints(cellvalues)
            dx = getdetJdV(cellvalues, qp)
            coords_qp = spatial_coordinate(cellvalues, qp, coords)
            s = coord_transform(coords_qp[2])
            coef[1, 1] = s
            coef[2, 2] = 1/s
            for i in 1:ndofs_c
                v = shape_value(cellvalues, qp, i)
                grad_v = shape_gradient(cellvalues, qp, i)
                for j in 1:ndofs_c
                    u = shape_value(cellvalues, qp, j) 
                    grad_u = shape_gradient(cellvalues, qp, j)
                    Ae[i, j] += (dot(coef*grad_u, grad_v) - (2*im*qm*s*grad_u[1]*v) + (qm^2*s*u*v)) * dx
                end
            end
        end
        
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A
end

function assemble_b(cellvalues::CellValues, B::SparseMatrixCSC, dofhandler::DofHandler)
    ndofs_c = ndofs_per_cell(dofhandler)
    assembler = start_assemble(B)
    
    # Preallocate the local stiffness matrix
    Be = zeros(ComplexF64, ndofs_c, ndofs_c)

    for cell in CellIterator(dofhandler)
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        fill!(Be, 0)
        coords = getcoordinates(cell)
        
        for qp in 1:getnquadpoints(cellvalues)
            dx = getdetJdV(cellvalues, qp)
            coords_qp = spatial_coordinate(cellvalues, qp, coords)
            coef = coord_transform(coords_qp[2])
            for i in 1:ndofs_c
                v = shape_value(cellvalues, qp, i)
                for j in 1:ndofs_c
                    u = shape_value(cellvalues, qp, j)
                    Be[i, j] += (coef * u * v) * dx
                end
            end
        end
        assemble!(assembler, celldofs(cell), Be)
    end

    return B
end

function main()
    d = 2*pi
    n = 20
    b = 2
    δ = 0.5
    grid = setup_grid(d, n, b+δ)

    # Define the interpolation: linear lagrange
    ip = Lagrange{RefTriangle, 1}()

    cv = setup_fevs(ip)
    dh = setup_dofs(grid, ip)
    ch = setup_bdcs(dh, d)

    # Allocate the matrices
    A = allocate_matrices(dh, ch)
    B = allocate_matrices(dh, ch)

    # Brillouin zone α
    α = collect(range(-pi/d, pi/d, 80))

    # Solve the eigenvalue problem
    k = get_dispersion(cv, dh, ch, A, B, α, 5)

    # Plot the dispersion curves
    p1 = plot(α, real(k), title = "Real part")
    p2 = plot(α, imag(k), title = "Imaginary part")
    plot(p1, p2, layout = (1, 2), legend = false)
end

main()