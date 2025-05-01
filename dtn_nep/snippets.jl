#----------------------------------------------------------
# Finite element discretization
#----------------------------------------------------------

"""
    assemble_A0(cv::CellValues, dh::DofHandler, A₀::SparseMatrixCSC, medium::Function, k)

Assemble the zero order term ``\\mathbf{A}_{0}`` in the Nonlinear eigenvalue problem.

``(\\mathbf{A}_{0} + \\alpha \\mathbf{A}_{1} + \\alpha^2 \\mathbf{A}_{2} + \\mathbf{F}(\\alpha)) \\mathbf{u} = 0.``

# Argument

- `cv`: CellValues (a notion in Ferrite.jl)
- `dh`: DofHandler (a notion in Ferrite.jl)
- `A₀`: an empty sparse pattern for A₀
- `medium`: refractive index function which describes the properties of the medium
- `k`: positive wavenumber
"""
function assemble_A0(cv::CellValues, dh::DofHandler, A₀::SparseMatrixCSC, medium::Function, k)
    # Allocate the local matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)
    
    # Create an assembler A₀ 
    assembler = start_assemble(A₀)
    
    # Loop over all cells
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)

        # Reset local stiffness matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)
        coords = getcoordinates(cell)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            n = medium(coords_qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                ∇v = shape_gradient(cv, qp, i)

                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble local stiffness matrix
                    Ae[i, j] += (∇u[1] * ∇v[1] + ∇u[2] * ∇v[1] - (k^2) * n * u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₀ 
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₀
end

"""
    assemble_A1(cv::CellValues, dh::DofHandler, A₁::SparseMatrixCSC)

Assemble the first order term ``\\mathbf{A}_{1}`` in the Nonlinear eigenvalue problem.

`` (\\mathbf{A}_{0} + \\alpha \\mathbf{A}_{1} + \\alpha^2 \\mathbf{A}_{2} + \\mathbf{F}(\\alpha)) \\mathbf{u} = 0. ``
"""
function assemble_A1(cv::CellValues, dh::DofHandler, A₁::SparseMatrixCSC)
    # Allocate the local stiffness matrix 
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)

    # Create an assembler 
    assembler = start_assemble(A₁)
    
    # Loop over all cells
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        
        # Reset local stiffness matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)
        
        # Loop over quadrature points 
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs 
                v = shape_value(cv, qp, i)

                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble local stiffness matrix 
                    Ae[i, j] += (-2im * ∇u[1] * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₁
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₁
end

"""
    assemble_A2(cv::CellValues, dh::DofHandler, A₂::SparseMatrixCSC)

Assemble the second order term in the Nonlinear eigenvalue problem.

`` (\\mathbf{A}_{0} + \\alpha \\mathbf{A}_{1} + \\alpha^2 \\mathbf{A}_{2} + \\mathbf{F}(\\alpha)) \\mathbf{u} = 0. ``
"""
function assemble_A2(cv::CellValues, dh::DofHandler, A₂::SparseMatrixCSC)
    # Allocate the local matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)

    # Create an assembler A₂
    assembler = start_assemble(A₂)

    # Loop over all cells 
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell 
        reinit!(cv, cell)
        
        # Reset local stiffness matrix to 0.0 + 0.0im 
        fill!(Ae, 0.0 + 0.0im)
        
        # Loop over quadrature points 
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                
                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    
                    # Assemble local stiffness matrix 
                    Ae[i, j] += (u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₂ 
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₂
end 

"""
    assemble_tbc(fv::FacetValues, dh::DofHandler, F::SparseMatrixCSC, facetset, dofsDtN, N, k, α)

Assemble the TBC matrix for the nonlinear eigenvalue problem.

# Arguments

- `fv`: `FacetValues`
- `dh`: `DofHandler`
- `F`: an empty sparse pattern for the TBC matrix
- `facetset`: `Facets` on the artificial boundary
- `dofsDtN`: Degrees of Freedom associated to the artificial boundary
- `N`: the number of the terms in the truncated DtN map
- `k`: the wavenumber
- `α`: the quasimomentum
"""
function assemble_tbc(fv::FacetValues, dh::DofHandler, F::SparseMatrixCSC, facetset, dofsDtN, N, k, α)
    # Allocate the vector Θ
    Θ = sparsevec(dofsDtN, zeros(ComplexF64, length(dofsDtN)), ndofs(dh))
    
    # Loop over truncated terms
    for n in -N:N 
        # Reset the vector Θ to zero
        fill!(Θ, 0.0 + 0.0im)
        
        # Compute βₙ 
        βₙ = beta_n(k, α, n)

        # Compute the vector Θ (Fourier coefficients and their conjugates)
        compute_coef!(fv, dh, facetset, Θ, n)
        
        # Assemble the TBC matrix
        for i in Θ.nzind, j in Θ.nzind
            v = im * βₙ * Θ[i] * conj(Θ[j])/(2π)
            Ferrite.addindex!(F, v, i, j)
        end
    end
    
    return F
end

function beta_n(k, α, n)
    αₙ = α + complex(n)
    βₙ = sqrt(complex(k)^2 - αₙ^2)
    
    return βₙ
end

#----------------------------------------------------------
# Nonlinear eigenvalue problem
#----------------------------------------------------------
"""
    Nep

Nep is used to store the matrices assocaited to the quadratic part.

    (A₀ + α A₁ + α² A₂ + F(α))u = 0 
"""
struct Nep 
    A₀::SparseMatrixCSC
    A₁::SparseMatrixCSC
    A₂::SparseMatrixCSC
end

# TODO: construct nonlinear eigenvalue problem u::Nep by u = Nep(A0, A1, A2)
# TODO: function (L::Nep)(\alpha::Complex)

"""
    (nep::Nep)(α)

TBW
"""
function (nep::Nep)(α)
    # TODO:Generate TBC matrix F(α)

    N = nep.A₀ + α * nep.A₁ + (α^2) * nep.A₂
    return N
end

"""
    Nnep{T}

Contain the information of the nonlinear eigenvalue problem (NEP) generated by 
the DtN-FEM. By the finite element discretization, we obtain

`` (\\mathbf{A}_{0} + \\alpha \\mathbf{A}_{1} + \\alpha^2 \\mathbf{A}_{2} + \\mathbf{F}(\\alpha)) \\mathbf{u} = 0. ``

We split the NEP into two parts: quadratic part and fully nonlinear part. For 
the quadraic part, we can decouple the `α` from the NEP easily. So the struct 
`Nnep` stores these three matrices directly. For the fully nonlinear part, `Nnep` 
just stores neccessary information for its finite element discretization.

# Fields

- `A₀`, `A₁`, `A₂`: the matrices associated to the quadratic part
- `fv`: `FacetValues` used for the integral on boundaries in Ferrite.jl
- `dh`: `DofHandler`
- `cst`: `ConstraintHandler` used for the information on boundary conditions
- `facetset`: the set of `Facets` on the artificial boundary
- `dofsDtN`: Degrees of Freedom associated to the artificial boundary
- `N`: the number of the terms in the truncated DtN map
- `k`: the wavenumber
"""
struct Nnep{T}
    A₀::SparseMatrixCSC{T, Int}
    A₁::SparseMatrixCSC{T, Int}
    A₂::SparseMatrixCSC{T, Int}
    fv::FacetValues
    dh::DofHandler
    cst::ConstraintHandler
    facetset
    dofsDtN
    N::Int64
    k
end

function (nep::Nnep{T})(α::T) where T
    # Allocate the TBC matrix
    F = allocate_stiff_matrix(nep.dh, nep.cst, nep.dofsDtN)

    # TODO:Assemble the TBC matrix 
    F = assemble_tbc(nep.fv, nep.dh, F, nep.facetset, nep.dofsDtN, nep.N, nep.k, α)
    
    # Impose the boundary conditions
    apply!(F, nep.cst)

    return nep.A₀ + α * nep.A₁ + (α^2) * nep.A₂ + F
    
end

function cim(ctr::Cim.AbstractContour, nep::Nnep, d::Int, l::Int; n=50, tol=1e-12)
    # Input validation
    d > 0 || throw(ArgumentError("d must be positive"))
    l > 0 || throw(ArgumentError("l must be positive"))

    # Get the quadrature points
    pts = get_quadpts(ctr, n)

    # Preallocate arrays
    A0 = zeros(ComplexF64, d, l)
    A1 = zeros(ComplexF64, d, l)
    Vhat = randn(ComplexF64, d, l)

    # Compute A0 and A1 with trapezoid rule
    for j in 1:pts.N
        z = complex(pts.nodes[j, 1], pts.nodes[j, 2])
        z_prime = complex(pts.nodes_prime[j, 1], pts.nodes_prime[j, 2])
        invNEP_Vhat = nep(z) \ Vhat
        A0 .+= invNEP_Vhat * z_prime
        A1 .+= invNEP_Vhat * z * z_prime
    end
    A0 ./= (pts.N*im)
    A1 ./= (pts.N*im)

    # Compute the SVD of A0
    (V, Sigma, W) = svd(A0)

    # Handle rank deficiency
    if isempty(Sigma)
        @warn "No eigenvalues found!"
        return ComplexF64[]
    end

    # Determine the number of nonzero singular values 
    k = count(Sigma ./ Sigma[1] .> tol)

    # Compute the matrix B 
    Vk = V[:,1:k]
    Sigk = Sigma[1:k]
    Wk = W[:,1:k]

    # Diagonal is more efficient
    B = (Vk' * A1 * Wk) * Diagonal(1 ./ Sigk)

    # Compute the eigenvalues of B 
    lambda = eigvals(B)

    # Avoid spurious eigenvalues
    filter!(λ -> Cim.is_inside(λ, ctr), lambda)

    return lambda
end
