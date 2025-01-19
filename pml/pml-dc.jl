#=
# pml-dc.jl
#
# Use PML to truncate the unbounded domain.
# 
# Author: Jiayi Zhang
# Date:   17/01/2025
# TODO: 1. periodic mesh in setup_grid()
#       2. generate A and B
=#

using Ferrite, FerriteGmsh
using SparseArrays, StaticArrays, Arpack
using Gmsh

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
    setup_grid(d, n, b, δ, lc, lp)

Generate the mesh by julia interface for Gmsh and read the `.msh` file by `FerriteGmsh`

# Arguments

- `d`: the period of the periodic grating
- `n`: the number of points on the profile of the grating
- `b`: the height of th truncated domain
- `lc`: target mesh size close to the point

- `δ`: the width of the pml layer
- `lp`: target mesh size in the pml layer

# Points

```
```
and
```
n+4 ------------- n+3
 |                 |
n+2 ------------- n+1
 |                 |
 1  -------------  n
```

# Lines

```
```
and
```
 .  ---- n+4 ----  .
n+5               n+3
 .  ---- n+1 ----  .
n+2                n 
 .  -- 1 ~ n-1 --  .
```

# Examples

```julia
d=2*pi, n=10, b=2
setup_grid(2*pi, 10, 2)
```

```julia
d=2*pi, n=10, b=2, δ=0.5
setup_grid(2*pi, 20, 2, 0.5)
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
    points = @MVector zeros(Int32, n+2)

    points[1] = gmsh.model.geo.addPoint(0, 0, 0, lc)
    points[n] = gmsh.model.geo.addPoint(d, 0, 0, lc)

    for i in 2:n-1
        points[i] = gmsh.model.geo.addPoint((i-1)*h, sin((i-1)*h), 0, lc) 
    end

    points[end-1] = gmsh.model.geo.addPoint(d, b, 0, lc) 
    points[end] = gmsh.model.geo.addPoint(0, b, 0, lc)

    # Add the lines
    lines = @MVector zeros(Int32, n+5)

    for i in 1:n+1
        lines[i] = gmsh.model.geo.addLine(points[i], points[i+1])
    end

    lines[end] = gmsh.model.geo.addLine(points[end], points[1])

    # Add the loop
    loop = gmsh.model.geo.addCurveLoop(lines)

    # Add the surface
    surf = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.synchronize()

    # Create the Physical Groups
    gmsh.model.addPhysicalGroup(1, [lines[end-1]], -1, "Top")
    gmsh.model.addPhysicalGroup(1, [lines[end]], -1, "Left")
    gmsh.model.addPhysicalGroup(1, [lines[end-2]], -1, "Right")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Omega")

    # Set periodic mesh
    transform = @SVector [1, 0, 0, d, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [rb], [lb], transform)

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

function setup_grid(d, n, b, δ, lc=0.5, lp=0.1) 
    h = d/(n-1)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    points = @MVector zeros(Int32, n+4)
    
    points[1] = gmsh.model.geo.addPoint(0, 0, 0, lc)
    points[n] = gmsh.model.geo.addPoint(d, 0, 0, lc)
    
    for i in 2:n-1
        points[i] = gmsh.model.geo.addPoint((i-1)*h, sin((i-1)*h), 0, lc) 
    end

    points[n+1] = gmsh.model.geo.addPoint(d, b, 0, lp)
    points[n+2] = gmsh.model.geo.addPoint(0, b, 0, lp)
    points[n+3] = gmsh.model.geo.addPoint(d, b+δ, 0, lp)
    points[n+4] = gmsh.model.geo.addPoint(0, b+δ, 0, lp)

    # Add the lines
    lines = @MVector zeros(Int32, n+5)

    for i in 1:n+1
        lines[i] = gmsh.model.geo.addLine(points[i], points[i+1])
    end

    lines[n+2] = gmsh.model.geo.addLine(points[n+2], points[1])
    lines[n+3] = gmsh.model.geo.addLine(points[n+1], points[n+3])
    lines[n+4] = gmsh.model.geo.addLine(points[n+3], points[n+4])
    lines[n+5] = gmsh.model.geo.addLine(points[n+4], points[n+2])

    lb = gmsh.model.geo.addLine(points[n+4], points[1])
    rb = gmsh.model.geo.addLine(points[n], points[n+3])

    # Add the loop
    lp1 = gmsh.model.geo.addCurveLoop(lines[1:n+2])
    lp2 = gmsh.model.geo.addCurveLoop([-lines[n+1], lines[n+3], lines[n+4], lines[n+5]])

    # Add the surface
    domain = gmsh.model.geo.addPlaneSurface([lp1])
    pml = gmsh.model.geo.addPlaneSurface([lp2])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the Physical Groups
    gmsh.model.addPhysicalGroup(1, [lines[n+4]], -1, "Top")
    gmsh.model.addPhysicalGroup(1, [lb], -1, "Left")
    gmsh.model.addPhysicalGroup(1, [rb], -1, "Right")
    gmsh.model.addPhysicalGroup(2, [domain, pml], -1, "Omega")

    # Set periodic mesh
    transform = @SVector [1, 0, 0, d, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [rb], [lb], transform)

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

function setup_bdcs(dofhandler::DofHandler)
    csthandler = ConstraintHandler(dofhandler)

    # Set periodic boundary condition
    periodic_faces = collect_periodic_facets(dofhandler.grid, "Right", "Left", x -> x + Vec{2}((2*pi, 0.0)))
    pbc = PeriodicDirichlet(:u, periodic_faces)
    add!(csthandler, pbc)

    # Set Dirichlet boundary condition on the "Top"
    dbc = Dirichlet(:u, getfacetset(dofhandler.grid, "Top"), x -> 0)
    add!(csthandler, dbc)

    close!(csthandler)
    return csthandler
end

function allocate_matries(dofhandler::DofHandler, csthandler::ConstraintHandler)
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
    dispersion = @MMatrix zeros(m, neigs)

    # Assemble the matrix B
    B = assemble_b(cellvalues, B, dofhandler)
    apply!(B, csthandler)

    # For fixed α
    for (i, qm) in enumerate(qms)
        # Assemble the matrix A(α)
        A = assemble_a(cellvalues, A, dofhandler, qm)
        apply!(A, csthandler)

        # solve the general eigenvalue problem
        λ, _ = eigs(A, B, nev = neigs, which = :LM)
        disperion[i, :] = λ
        fill!(A, 0.0)
    end
    
    return dispersion
end

function assemble_a(cellvalues::CellValues, A::SparseMatrixCSC, dofhandler::DofHandler, qm)
    ndofs_c = ndofs_per_cell(dofhandler)
    assembler = start_assemble(A)
    coef = @MMatrix zeros(2, 2)
    
    # Preallocate the local stiffness matrix
    Ae = zeros(ComplexF64, ndofs_c, ndofs_c)

    for cell in CellIterator(dofhandler)
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        fill!(Ae, 0)

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
    ch = setup_bdcs(dh)

    # Allocate the matrices
    A = allocate_matrices(dh, ch)
    B = allocate_matrices(dh, ch)

    # Brillouin zone α
    α = collect(range(-pi/d, pi/d, 10))

    #
    k = get_dispersion(cv, dh, ch, A, B, α, 5)

    # Plot the dispersion curves
    @show k

    return nothing
end

main()