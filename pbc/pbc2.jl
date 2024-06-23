"""
    pbc2.jl -- This file works on v0.3.14
    -------------------------------------------------------
    PROBLEMS:
                 Gammab
           ----------------------
           |                    |
    Gamma2 |                    | Gamma1
           |                    |
           ----------------------
                 Gamma

    Δu = -2cos(x1)sin(x2)
    u = 0 on Gamma and Gammab
    periodic boundary condition on Gamma1 and Gamma2

    The analytic solution of u is u(x1, x2) = cos(x1)sin(x2)

    GOALS:
          1. generate grid by using Gmsh
          2. implement the periodic boundary condition

    REFERENCES:
          1. Heat equation in Ferrite document: https://ferrite-fem.github.io/Ferrite.jl/stable/examples/heat_equation/
          2. Helmholtz equation in Ferrite document: https://ferrite-fem.github.io/Ferrite.jl/stable/examples/helmholtz/
          3. Stokes flow: https://ferrite-fem.github.io/Ferrite.jl/stable/examples/stokes-flow/
"""
# If you use gradient(), then you should add using Tensor
using Ferrite, SparseArrays, FerriteGmsh, Gmsh

# R.H.S.
function rhs(x::Vec{2, T}) where {T}
    return -2*cos(x[1])*sin(x[2])
end

# Generate the mesh by Gmsh
function setup_grid(h=0.05)
    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, h)
    p2 = gmsh.model.geo.addPoint(2*pi, 0, 0, h)
    p3 = gmsh.model.geo.addPoint(2*pi, pi, 0, h)
    p4 = gmsh.model.geo.addPoint(0, pi, 0, h)
  
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
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Omega")
    
    # Set Periodic boundary condition
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], [1, 0, 0, 2*pi, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)  

    # Save the mesh, and read the .msh file in Ferrite
    path = joinpath(pwd(), "mesh.msh") 
    gmsh.write(path)
    grid = togrid(path)

    # Finalize the Gmsh library 
    gmsh.finalize()

    return grid
end

#
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

function main()
    # Setup the grid
    h = 0.5 
    grid = setup_grid(h)

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
    periodic_faces = collect_periodic_faces(grid, "Gamma1", "Gamma2", x -> x + Vec{2}((2*pi, 0.0)))
    pbc = PeriodicDirichlet(:u, periodic_faces)
    add!(ch, pbc)
    # Dirichlet BC
    dbc = Dirichlet(:u, union(getfaceset(grid, "Gamma"), getfaceset(grid, "Gammab")), x -> 0)
    add!(ch, dbc)
    close!(ch)
    
    K = create_sparsity_pattern(dh, ch)
    K, f = doassemble(cellvalues, K, dh)
    apply!(K, f, ch)
    u = K \ f;
    apply!(u, ch)
    
    # Write to VTK
    vtk_grid("pbc2", dh) do vtk 
        vtk_point_data(vtk, dh, u)
    end
end

main()
