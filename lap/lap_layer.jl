using Ferrite
using Nmrc
using SparseArrays

# functions in Nmrc we need in this script: 
#     periodic_cell, setup_vals, setup_dh, setup_bcs,
#     allocate_stiff_matrix, assemble_A0, assemble_A1, assemble_A2, 
#     compute_coef!, dofs_on_dtn, Incident, assemble_load, assemble_tbc

# For CIM
using Random, LinearAlgebra, Cim

# ---------------------------------------------------------
# Define functions
# ---------------------------------------------------------

# The refractive index

function medium(x)
    if x[2] < 1.0
        return 3.9
    else
        return 1.0
    end
end

# TODO: We have defined csqrt_negimag in Nmrc.
#       We also need a beta_n for complex α
function my_sqrt(z::Complex{T}) where T<:AbstractFloat
    # Handle special case of zero input
    if z == zero(z)
        return zero(z)
    end
    
    # Rotate by -π/2 (multiply by -i)
    rot_z = Complex(imag(z), -real(z))
    
    # Compute standard square root
    sqrt_rot = sqrt(rot_z)
    
    # Rotate by π/4 (multiply by exp(im * π/4))
    e = 1/sqrt(2)
    r = (real(sqrt_rot) - imag(sqrt_rot)) * e
    i = (real(sqrt_rot) + imag(sqrt_rot)) * e
    result = Complex(r, i)
    
    return result
end

# Convenience method for other number types
my_sqrt(z::Complex) = my_sqrt(float(z))

function beta_n(k, α, n)
    αₙ = α + complex(n)
    βₙ = my_sqrt(complex(k)^2 - αₙ^2)
    
    return βₙ
end

# Assemble TBC part for nonlinear eigenvalue problem
function assemble_tbc(fv::FacetValues, dh::DofHandler, F::SparseMatrixCSC, facetset, dofsDtN, N, k, α)
    # Allocate the vector Θ
    Θ = sparsevec(dofsDtN, zeros(ComplexF64, length(dofsDtN)), ndofs(dh))
    
    # Loop over truncated terms
    for n in -N:N 
        # Reset the vector Θ to zero
        fill!(Θ, zero(ComplexF64))
        
        # Compute βₙ
        βₙ = beta_n(k, α, n)
        
        # Compute the vector Θ (Fourier coefficients and their conjugates)
        Nmrc.compute_coef!(fv, dh, facetset, Θ, n)
        
        # Assemble the TBC matrix
        for i in Θ.nzind, j in Θ.nzind
            v = im * βₙ * Θ[i] * conj(Θ[j]) / (2π)
            Ferrite.addindex!(F, v, i, j)
        end
    end
    
    return F
end

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
    
    # Assemble the TBC matrix
    F = assemble_tbc(nep.fv, nep.dh, F, nep.facetset, nep.dofsDtN, nep.N, nep.k, α)
    
    # Impose the boundary conditions
    apply!(F, nep.cst)
    
    return nep.A₀ + α * nep.A₁ + (α^2) * nep.A₂ + F
end

function new_cim(ctr::Cim.AbstractContour, nep::Nnep, d::Int, l::Int; n=50, tol=1e-12)
    # Input validation
    d > 0 || throw(ArgumentError("d must be positive"))
    l > 0 || throw(ArgumentError("l must be positive"))

    # Get the quadrature points
    pts = get_quadpts(ctr, n)

    # Preallocate arrays
    A0 = zeros(ComplexF64, d, l)
    A1 = zeros(ComplexF64, d, l)
    Random.seed!(10);
    Vhat = randn(ComplexF64, d, l)

    # Compute A0 and A1 with trapezoid rule
    for j in 1:pts.N - 1
        z = complex(pts.nodes[j, 1], pts.nodes[j, 2])
        z_prime = complex(pts.nodes_prime[j, 1], pts.nodes_prime[j, 2])
        invNEP_Vhat = nep(z) \ Vhat
        A0 .+= invNEP_Vhat * z_prime
        A1 .+= invNEP_Vhat * z * z_prime
    end
    A0 ./= (im * (pts.N - 1))
    A1 ./= (im * (pts.N - 1))

    # Compute the SVD of A0
    (V, Sigma, W) = svd(A0)
    @show Sigma

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

    # Compute the eigenvalues λ and the corresponding eigenvectors s
    λ, s = eigen(B)

    # Check if the eigenvalue lies inside the contour
    # filter!(λ -> Cim.is_inside(λ, ctr), lambda)
    
    # Store the indices of eigenvalues inside the contour
    good = Int[]
    
    for i in eachindex(λ)
        if Cim.is_inside(λ[i], ctr) == true
            push!(good, i)
        end
    end
    

    λ = @view λ[good]
    u = @view s[:, good]
    u = Vk * u

    return λ, u
end

# Assemble the stiffness matrix for layer problems
function assemble(cv::CellValues, dh::DofHandler, A::SparseMatrixCSC, n::Function, inc::Incident)
    k = get_wavenumber(inc)
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
        
        # Reset local stiffness matrix to 0.0 + 0.0
        fill!(Ae, 0.0 + 0.0im)
        
        # Get the coordinates of this cell
        coords = getcoordinates(cell)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Get the coordinates of the quadrature point and 
            # evaluate the refractive index at the quadrature point
            coords_qp = spatial_coordinate(cv, qp, coords)
            ri = n(coords_qp)
            
            # Loop over test shape functions
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                ∇v = shape_gradient(cv, qp, i)

                # Loop over trial shape functions
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble the local stiffness matrix
                    Ae[i, j] += (∇u ⋅ ∇v - 2im * α * ∇u[1] * v - ((k^2) * ri - (α^2)) * u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A
end

# Compute the integral 
function integral_uv̄(cv::CellValues, dh::DofHandler, u, v)
    # Check if the elements in u and v have the same type
    (eltype(u) == eltype(v)) || throw(ArgumentError("elements in u and v should have the same type!")) 
    
    rst = zero(eltype(u))
    
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        
        # Get the corresponding cell id of the current cell cache
        ci = cellid(cell)
        
        # Get the DoFs of the cell `ci`
        cd = celldofs(dh, ci)

        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            u_qp = function_value(cv, qp, u, cd)
            v_qp = function_value(cv, qp, v, cd)
            rst += u_qp * conj(v_qp) * dx
        end
    end
    
    return rst
end

# Compute the integral
function integral_nuv̄(cv::CellValues, dh::DofHandler, n::Function, u, v)
    # Check if the elements in u and v have the same type
    (eltype(u) == eltype(v)) || throw(ArgumentError("elements in u and v should have the same type!")) 
    
    rst = zero(eltype(u))
    
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        
        # Get the corresponding cell id of the current cell cache
        ci = cellid(cell)
        
        # Get the DoFs of the cell `ci`
        cd = celldofs(dh, ci)
        
        # Get the coordinates of this cell
        coords = getcoordinates(cell)

        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Get the coordinate of the quadrature point 
            coords_qp = spatial_coordinate(cv, qp, coords)

            # Evaluate the function value of n
            ri = n(coords_qp)

            u_qp = function_value(cv, qp, u, cd)
            v_qp = function_value(cv, qp, v, cd)
            rst += ri * u_qp * conj(v_qp) * dx
        end
    end
    
    return rst
end

# Compute the integral
function integral_∂₁uv̄(cv::CellValues, dh::DofHandler, u, v)
    # Check if the elements in u and v have the same type
    (eltype(u) == eltype(v)) || throw(ArgumentError("elements in u and v should have the same type!"))  

    rst = zero(eltype(u))
    
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        
        # Get the corresponding cell id of the current cell cache
        ci = cellid(cell)
        # Get the DoFs of the cell `ci`
        cd = celldofs(dh, ci)

        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            ∂₁u_qp = function_gradient(cv, qp, u, cd)[1]
            v_qp = function_value(cv, qp, v, cd)
            rst += ∂₁u_qp * conj(v_qp) * dx
        end
    end
    
    return rst
end

# ---------------------------------------------------------
# Define constants and common variables
# ---------------------------------------------------------

# Wavenumber
k = 4.1

# The number of the truncated terms in Rayleigh expansions
N = 20

# Mesh for the periodic cell
grid = periodic_cell(lc = 0.01, period = 2π, height = 2.5)

# Interpolation for FEM
ip = Lagrange{RefTriangle, 1}()

# Set up CellValues and FacetValues
cv, fv = setup_vals(ip)

# Set up the DofHandler
dh = setup_dh(grid, ip)

# Set up the ConstraintHandler: 
# Dirichlet boundary conditon on the bottom and periodic conditon on the left and right
cst = setup_bcs(dh)

# Extract the facetset on the "top" boundary
top = getfacetset(grid, "top")

# Extract the DoFs on the "top" boundary
dofsDtN = dofs_on_dtn(dh, :u, top)

# ---------------------------------------------------------
# Step 1: Get the exceptional value and corresponding 
#         eigenfunctions
# ---------------------------------------------------------

# Allocate A₀, A₁, A₂ for the quadratic part
A₀ = allocate_stiff_matrix(dh, cst, dofsDtN)
A₁ = allocate_stiff_matrix(dh, cst, dofsDtN)
A₂ = allocate_stiff_matrix(dh, cst, dofsDtN)

# Allocate A₀, A₁, A₂
A₀ = assemble_A0(cv, dh, A₀, medium, k)
A₁ = assemble_A1(cv, dh, A₁)
A₂ = assemble_A2(cv, dh, A₂)

# Impose the constraint on A₀, A₁, and A₂
apply!(A₀, cst)
apply!(A₁, cst)
apply!(A₂, cst)

# The size of the matrix needed by CIM
d = size(A₀, 1)

# Construct the nonlinear eigenvalue problem
L = Nnep(A₀, A₁, A₂, fv, dh, cst, top, dofsDtN, N, k); 

# Define the contour around the point obtained by PML result
cir = Cim.circle([0.38, 0.0], 0.01) 

# Solve the nonlinear eigenvalue problem by CIM
λ, v = new_cim(cir, L, d, 5; n = 50, tol = 1e-7)

# Get the eigenvalues and the corresponding eigenfunctions
α = λ[1]
ϕ = v[:, 1]

# Apply the constraint on the eigenfunction
apply!(ϕ, cst)

# ---------------------------------------------------------
# Step 2: Get the particular solution w.r.t. the exceptional 
#         value we got
# ---------------------------------------------------------

# Compute the incident angle in radians based on the 
# exceptional value in Step 1
θ = asin(real(α) / k)

inc = Incident(k, θ)

# Allocate matrices
A = allocate_stiff_matrix(dh, cst, dofsDtN)
F = allocate_stiff_matrix(dh, cst, dofsDtN)
f = zeros(ComplexF64, ndofs(dh))

# Assemble the load vector f
f = assemble_load(fv, dh, top, f, inc, 2.5)

# Assemble the matrix A
A = assemble(cv, dh, A, medium, inc)

# Assemble the TBC matrix
F = Nmrc.assemble_tbc(fv, dh, inc, top, F, N, dofsDtN)

# Add the TBC matrix to A
A = sub_preserve_structure(A, F)

# Impose the boundary condition
apply!(A, f, cst)

# Solve the linear system
u₀ = A \ f

# Apply the boundary condition again
# Note: I forgot this step in the first version
apply!(u₀, cst)

# ---------------------------------------------------------
# Step 3: Compute the coefficient according to the 
#         constraint proposed by G. Hu and A. Kirsch. 
#
# Note:!! We note that the constraint in G. Hu and A. Kirsch's paper 
#         holds on quasi-periodic function spaces. But our numerical 
#         solutions are all periodic due to issues of implementation. 
#         So we need to reformulate the constraint for a periodic setting. 
# Note:!! For layer problems, we need a similar constraint.
# ---------------------------------------------------------

# Compute constants in the formulation
s = sin(θ)
s² = sin(θ)^2

# There are four integrals to compute

# 1. ∂ϕ ̄ϕ
c1 = integral_∂₁uv̄(cv, dh, ϕ, ϕ)

# 2. n |ϕ|^2
c2 = real(integral_nuv̄(cv, dh, medium, ϕ, ϕ))

# 3. |ϕ|^2
c3 = integral_uv̄(cv, dh, ϕ, ϕ)

# 4. ∂u₀ ̄ϕ
c4 = integral_∂₁uv̄(cv, dh, u₀, ϕ)

# 5. n u₀ ̄ϕ
c5 = integral_nuv̄(cv, dh, medium, u₀, ϕ)

# 6. u₀ ϕ̄
c6 = integral_uv̄(cv, dh, u₀, ϕ)

# Compute the coefficient
coef = (-s * c4 + im * k * c5 - im * k * s² * c6) / (s * c1 - im * k * c2 + im * k * s² * c3)

# ---------------------------------------------------------
# Step 4: Get the numerical solution
# ---------------------------------------------------------

u = u₀ + coef .* ϕ

VTKGridFile("modulus_u", grid) do vtk
    write_solution(vtk, dh, abs.(u))
end;

VTKGridFile("modulus_utheta", grid) do vtk
    write_solution(vtk, dh, abs.(u₀))
end;

VTKGridFile("modulus_eigf", grid) do vtk
    write_solution(vtk, dh, abs.(ϕ))
end;


# # From this part, we know that only around 3% elements in abs.(u) are greater than 3.0
# function is_greater3(x)
#     if x > 3.0
#         return true
#     else
#         return false
#     end
# end
# 
# ũ = filter(is_greater3, abs.(u))

# length(ũ) / length(u)
# 
# VTKGridFile("real_u", grid) do vtk
#   write_solution(vtk, dh, real.(u))
# end;

# VTKGridFile("imag_u", grid) do vtk
#    write_solution(vtk, dh, imag.(u))
# end;