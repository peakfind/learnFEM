#=
# dtn.jl
# 
# Provide a framework for the DtN-FEM 
# 
# Author: Jiayi Zhang
# Date: 23/01/2025
=#
using Gmsh, Ferrite, FerriteGmsh, SparseArrays

# The incident field
struct Incident
    k::Float64 # wave number
    θ::Float64 # incident angle
    α::Float64 # α = k*sin(θ)
    β::Float64 # β = k*cos(θ)
end

function Incident(k, θ)
    α = k*sin(θ)
    β = k*cos(θ)
    Incident(k, θ, α, β)
end

# Generate the mesh by Gmsh
function setup_grid(height, lc=0.05)
    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(2*pi, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(2*pi, height, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, height, 0, lc)

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
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], [1, 0, 0, 2*pi, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    # Save the mesh and read the .msh file in Ferrite
    path = joinpath(pwd(), "layer.msh")
    gmsh.write(path)
    grid = togrid(path)

    # Finalize the Gmsh library
    gmsh.finalize()

    return grid
end

function setup_vals(interpolation)
    qr = QuadratureRule{RefTriangle}(2)
    qr_facet = FacetQuadratureRule{RefTriangle}(2)
    cellvalues = CellValues(qr, interpolation)
    facetvalues = FacetValues(qr_facet, interpolation)
    return cellvalues, facetvalues
end

function setup_dofs(grid::Grid, interpolation)
    dh = DofHandler(grid)
    add!(dh, :u, interpolation)
    close!(dh)
    return dh
end

function setup_bcs(dofhandler::DofHandler) 
    ch = ConstraintHandler(dofhandler)
    # Periodic boundary condition
    periodic_faces = collect_periodic_facets(dofhandler.grid, "Gamma1", "Gamma2", x -> x + Vec{2}((2*pi, 0.0)))
    pbc = PeriodicDirichlet(:u, periodic_faces)
    add!(ch, pbc)
    # Dirichlet boundary condition
    # dirichlet_faces = union(getfacetset(dofhandler.grid, "Gamma"), getfacetset(dofhandler.grid, "Gammab"))
    dbc = Dirichlet(:u, getfacetset(dofhandler.grid, "Gamma"), x -> 0)
    add!(ch, dbc)
    close!(ch)
    return ch
end

# TODO: Maybe find a more elegant way to extract the dofs on DtN
function dofs_on_dtn(dofhandler::DofHandler, field::Symbol, facetset)
    dtn_ch = ConstraintHandler(dofhandler)
    dbc = Dirichlet(field, facetset, x -> 0)
    add!(dtn_ch, dbc)
    close!(dtn_ch)
    return dtn_ch.prescribed_dofs
end

function allocate_stiff_matrix(dofhandler::DofHandler, csthandler::ConstraintHandler, dofs)
    sp = init_sparsity_pattern(dofhandler)
    add_cell_entries!(sp, dofhandler)
    # Use add_entry! for DtN term
    for i in dofs, j in dofs
        if abs(i - j) > 1
            Ferrite.add_entry!(sp, i, j)
        end
    end
    add_constraint_entries!(sp, csthandler)
    K = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return K
end

function doassemble(cellvalues::CellValues, facetvalues::FacetValues, K::SparseMatrixCSC, f, dh::DofHandler, alpha, beta, b)
    ndofs_c = ndofs_per_cell(dh) 
    assembler = start_assemble(K, f)
    
    # Preallocate the local stiffness matrix and local load vector
    Ke = zeros(ComplexF64, ndofs_c, ndofs_c)
    fe = zeros(ComplexF64, ndofs_c)
    
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        fill!(Ke, 0)
        # coords = getcoordinates(cell)

        for qp in 1:getnquadpoints(cellvalues)
            dx = getdetJdV(cellvalues, qp)
            # coords_qp = spatial_coordinate(cellvalues, qp, coords)
            # ri = q(coords_qp)
            for i in 1:ndofs_c
                v = shape_value(cellvalues, qp, i)
                grad_v = shape_gradient(cellvalues, qp, i)
                for j in 1:ndofs_c
                    u = shape_value(cellvalues, qp, j)
                    grad_u = shape_gradient(cellvalues, qp, j)
                    Ke[i, j] += (dot(grad_u, grad_v) - (im*2*alpha*grad_u[1]*v) - (beta^2)*u*v) * dx  
                end
            end
        end
        assemble!(assembler, celldofs(cell), Ke)
    end

    for facet in FacetIterator(dh, getfacetset(dh.grid, "Gammab"))
        # Reinitialize facetvalues for this cell
        reinit!(facetvalues, facet)
        fill!(fe, 0)
        
        for qp in 1:getnquadpoints(facetvalues)
            ds = getdetJdV(facetvalues, qp)
            for i in 1:getnbasefunctions(facetvalues)
                v = shape_value(facetvalues, qp, i)
                fe[i] += (-2*im*beta*exp(-im*beta*b))*v*ds
            end
        end
        assemble!(f, celldofs(facet), fe)
    end

    return K, f
end

function tbc_matrix(truc, dofs)
end

function main()
    # Setup the grid
    grid = setup_grid(5.0, 0.3)
    
    # Setup the parameters:k = 1.0, θ = pi/3, N = 15
    inc = Incident(1.0, pi/3)
    N = 15

    # Interpolation and quadrature rules
    ip = Lagrange{RefTriangle, 1}()

    # Set the FE values
    cv, fv = setup_vals(ip)

    # Set the DofHandlers
    dh = setup_dofs(grid, ip)

    # Set the boundary conditions
    ch = setup_bcs(dh)

    # Extract the dofs related with DtN term
    dofs_dtn = dofs_on_dtn(dh, :u, getfacetset(grid, "Gammab"))

    # Allocate the stiffness matrix and the load vector
    K = allocate_stiff_matrix(dh, ch, dofs_dtn)
    f = zeros(ComplexF64, ndofs(dh))

    # Assemble the stiffness matrix (do not contain the DtN term) and the load vector
    # Note that the DtN term should be computed on global level
    K,f = doassemble(cv, fv, K, f, dh, inc.α, inc.β, 10.0)

    # apply the DtN term
    # K_tbc = tbc_matrix()
    # K .-= K_tbc;

    # Apply the constraints
    apply!(K, f, ch)

    u = K \ f;
    apply!(u, ch)
    
    # Write to .vtk file
    VTKGridFile("real_u", dh) do vtk
        write_solution(vtk, dh, real(u))
    end
    VTKGridFile("imag_u", dh) do vtk
        write_solution(vtk, dh, u)
    end
end

main()
