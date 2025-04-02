# dtn.jl
# Provide a framework for the DtN-FEM 

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
    return Incident(k, θ, α, β)
end

function get_alpha(inc::Incident)
    return inc.α
end

function get_beta(inc::Incident)
    return inc.β
end

function beta_n(inc::Incident, n)
    αₙ = inc.α + n
    
    if k > abs(αₙ) 
        βₙ = Complex(sqrt(k^2 - αₙ^2))
    else
        βₙ = im * sqrt(αₙ^2 - k^2)
    end
    
    return βₙ
end

"""
    periodic_cell(lc=0.05; height)

Generate a mesh for the periodic cell with period 2π.
"""
function periodic_cell(lc=0.05; height)
    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(2π, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(2π, height, 0, lc)
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
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom") 
    gmsh.model.addPhysicalGroup(1, [l2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")

    # Set Periodic boundary condition
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], [1, 0, 0, 2π, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    # Save the mesh and read the .msh file in Ferrite
    grid = mktempdir() do dir 
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path)
    end

    # Finalize the Gmsh library
    gmsh.finalize()

    return grid
end

# setup cellvalues and facetvalues
function setup_vals(interpolation)
    qr = QuadratureRule{RefTriangle}(2)
    qr_facet = FacetQuadratureRule{RefTriangle}(2)
    cellvalues = CellValues(qr, interpolation)
    facetvalues = FacetValues(qr_facet, interpolation)
    return cellvalues, facetvalues
end

# Same as the setup_dof in fem/setups.jl
function setup_dofs(grid::Grid, interpolation)
    dh = DofHandler(grid)
    add!(dh, :u, interpolation)
    close!(dh)
    return dh
end

function setup_bcs(dofhandler::DofHandler) 
    ch = ConstraintHandler(dofhandler)

    # Periodic boundary condition
    periodic_faces = collect_periodic_facets(dofhandler.grid, "right", "left", x -> x + Vec{2}((2π, 0.0)))
    pbc = PeriodicDirichlet(:u, periodic_faces)
    add!(ch, pbc)

    # Dirichlet boundary condition
    dbc = Dirichlet(:u, getfacetset(dofhandler.grid, "bottom"), x -> 0)
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
        # if abs(i - j) > 1
            Ferrite.add_entry!(sp, i, j)
        # end
    end
    add_constraint_entries!(sp, csthandler)
    K = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return K
end

function assemble_A(cv::CellValues, dh::DofHandler, A::SparseMatrixCSC, inc::Incident)
    α = get_alpha(inc)
    # Allocate the local stiffness matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)
    
    # Create an assembler A
    assembler = start_assemble(A)
    
    # Loop over all cells 
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for current cell 
        reinit!(cv, cell)
        # Reset local stiffness matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)
        
        # Loop over quadrature points 
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                ∇v = shape_gradient(cv, qp, i)
                
                # Loop over trial shape functions
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble the local stiffness matrix
                    Ae[i, j] += (∇u ⋅ ∇v - 2im * α * ∇u[1] * v - (k^2 - α^2) * u * v) * dx
                end
            end
        end

        # Assemble Ae into A
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A
end

function assemble_load(fv::FacetValues, dh::DofHandler, facetset, f, inc::Incident, height)
    β = get_beta(inc)
    # Allocate the local load vector fe
    n_basefuncs = getnbasefunctions(fv)
    fe = zeros(ComplexF64, n_basefuncs)
    
    # Loop over all facets on the specific facetset
    for facet in FacetIterator(dh, facetset)
        # Update the fv to the correct facet
        reinit!(fv, facet)
        
        # Reset the local vector fe to zero
        fill!(fe, 0.0 + 0.0im)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(fv)
            ds = getdetJdV(fv, qp)
            
            # right hand side due to the incident field
            g = -2im * β * exp(-im * β * height)
            
            # Loop over test functions
            for i in 1:n_basefuncs
                v = shape_value(fv, qp, i)
                fe[i] += g * v * ds
            end
        end
        assemble!(f, celldofs(facet), fe)
    end
    
    return f
end

function assemble_tbc(fv::FacetValues, dh::DofHandler, inc::Incident, facetset, F, N, dofsDtN)
    # Allocate the vector Θ 
    Θ = sparsevec(dofsDtN, zeros(ComplexF64, length(dofsDtN)), ndofs(dh))
    
    # Loop over truncated terms
    for n in -N:N 
        # Reset the vector Θ to zero
        fill!(Θ, 0.0 + 0.0im)

        # Compute βₙ
        βₙ = beta_n(inc, n)

        # Compute the vector Θ (Fourier coefficients and its conjugate) 
        compute_coef!(fv, dh, facetset, Θ, n)
        
        # 
        for i in Θ.nzind, j in Θ.nzind
            v = im * βₙ * Θ[i] * conj(Θ[j])/(2π)
            Ferrite.addindex!(F, v, i, j)
        end
    end
    
    # F .*= im/(2π)
    
    return F
end

"""
    compute_coef!(fv::FacetValues, dh::DofHandler, facetset, Θ::SparseVector, n)

Compute Θⁿ on the `facetset`.
"""
function compute_coef!(fv::FacetValues, dh::DofHandler, facetset, Θ::SparseVector, n)
    # Allocate the local vector θ
    n_basefuncs = getnbasefunctions(fv)
    θ = zeros(ComplexF64, n_basefuncs)
    
    # Loop over all facets on the specific facetset
    for facet in FacetIterator(dh, facetset)
        # Update the fv to the correct facet
        reinit!(fv, facet)
        
        # Reset the local vector θ to zero
        fill!(θ, 0.0 + 0.0im)
        
        coords = getcoordinates(facet)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(fv)
            ds = getdetJdV(fv, qp)

            # Coordinate of the quadrature point
            coords_qp = spatial_coordinate(fv, qp, coords)
            
            # Modes: eⁱⁿˣ
            mode = exp(im * n * coords_qp[1])
            
            for i in 1:n_basefuncs
                ϕ = shape_value(fv, qp, i)
                θ[i] += ϕ * mode * ds
            end
        end
        assemble!(Θ, celldofs(facet), θ)
    end
    
    return Θ
end

function sub_preserve_structure(A::SparseMatrixCSC, B::SparseMatrixCSC)
    if A.colptr != B.colptr || A.rowval != B.rowval || size(A) != size(B)
        error("Matrices must have the same sparsity structure and dimensions.")
    end
    new_nzval = A.nzval - B.nzval
    return SparseMatrixCSC(size(A,1), size(A,2), A.colptr, A.rowval, new_nzval)
end