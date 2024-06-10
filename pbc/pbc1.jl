""" 
    pbc1.jl
    -------------------------------------------------------
    PROBLEMS 
                  Gammab
           --------------------
           |                  |
    Gamma2 |                  | Gamma1
           |                  |
           --------------------
                  Gamma
    
    Δu = -2cos(x1)sin(x2)
    u = 0 on Gamma and Gammab
    periodic boundary condition on Gamma1 and Gamma2

    The analytic solution of u is u(x1, x2) = cos(x1)sin(x2)
    
    GOALS:
          1. generate grid by using generate_grid() in Ferrite
          2. implement the periodic boundary condition 
    REFERENCE:
          1. Heat equation in Ferrite document: https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equations/
          2. Helmholtz equation in Ferrite document: https://ferrite-fem.github.io/Ferrite.jl/stable/examples/helmholtz/
          3. Periodic boundary condition: https://ferrite-fem.github.io/Ferrite.jl/stable/manual/boundary_conditions/#Periodic-boundary-conditions
"""

using Ferrite, SparseArrays

# R.H.S.
function rhs(x::Vec{2, T}) where {T} 
    return -2*cos(x[1])*sin(x[2])
end

# Generate grid
grid = generate_grid(Triangle, (20, 20), Tensor{1,2}((0.0,0.0)), Tensor{1,2}((2*pi,0.0)), Tensor{1,2}((2*pi,pi)), Tensor{1,2}((0.0,pi)))

# Create DofHandler
ip = Lagrange{2, RefTetrahedron, 1}()
qr = QuadratureRule{2, RefTetrahedron}(2)
cellvalues = CellScalarValues(qr, ip);

dh = DofHandler(grid)
add!(dh, :u, 1)
close!(dh);

# Create Constraint
ch = ConstraintHandler(dh);
# Periodic BC
periodic_faces = collect_periodic_faces(grid, "left", "right")
pbc = PeriodicDirichlet(:u, periodic_faces)
add!(ch, pbc)
# Dirichlet BC
dbc = Dirichlet(:u, union(getfaceset(grid, "top"), getfaceset(grid, "bottom")), x -> 0)
add!(ch, dbc)
close!(ch)

# Assemble
function doassemble(cellvalues::CellScalarValues, K::SparseMatrixCSC, dh::DofHandler)
    # Preallocate the global load vector of the linear system
    f = zeros(ndofs(dh))

    assembler = start_assemble(K, f)
    n_basefuncs = getnbasefunctions(cellvalues)

    # Preallocate the local stiffness matrix and local load vector
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    
    # Loop over all cells
    for cell in CellIterator(dh)
        fill!(Ke, 0)
        fill!(fe, 0)
        coords = getcoordinates(cell)

        reinit!(cellvalues, cell)

        #
        for q_point in 1:getnquadpoints(cellvalues)
            dx = getdetJdV(cellvalues, q_point)
            coords_qp = spatial_coordinate(cellvalues, q_point, coords)
            rhs_qp = rhs(coords_qp)
            
            for i in 1:n_basefuncs
                v = shape_value(cellvalues, q_point, i)
                grad_v = shape_gradient(cellvalues, q_point, i)
                fe[i] += (rhs_qp * v) * dx
                for j in 1:n_basefuncs
                    grad_u = shape_gradient(cellvalues, q_point, j)
                    Ke[i, j] += dot(grad_u, grad_v) * dx
                end
            end
        end
        
        assemble!(assembler, celldofs(cell), Ke, fe)
    end

    return K, f
end

# Solve the linear system
K = create_sparsity_pattern(dh, ch);
K, f = doassemble(cellvalues, K, dh)
apply!(K, f, ch)
u = K \ f;
apply!(u, ch)

# Write to VTK
vtk_grid("pbc1", dh) do vtk
    vtk_point_data(vtk, dh, u)
end
