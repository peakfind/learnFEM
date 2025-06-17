using Ferrite
# reuse some functions
using Nmrc: periodic_cell, setup_vals, setup_dh, setup_bcs, dofs_on_dtn, allocate_stiff_matrix,
            assemble_A0, assemble_A1, assemble_A2, compute_coef!
using Cim
using LinearAlgebra
using SparseArrays
using Random

include("snippets.jl")

# Step 1: Choose the wavenumber and the refractive index.
# Step 2: Generate the mesh.
# Step 3: Define the contour. We need to choose the contour carefully.

## Parameters

# Wavenumber for the example 1
k = 8.1

# Wavenumber for the example 2
# k = 4.1

# The nb of the truncated terms
N = 20

# Refractive index for the example 1
function medium(x)
    n = 1.0

    if x[2] < 1.0
        n = 3.0 + sin(2.0 * x[1])
    end
    
    return n
end

# Refractive index for the example 2
# function medium(x)
#     if x[2] < 1.0
#         return 3.9
#     else
#         return 1.0
#     end
# end

# Generate the mesh in a periodic cell
# grid = periodic_cell(lc=0.03, period=2π, height=2.0) # coarser mesh
grid = periodic_cell(lc=0.01, period=2π, height=2.0) # finer mesh

# Set up fevalues(CellValues and FacetValues), DofHandler, and ConstraintHandler
ip = Lagrange{RefTriangle, 1}()
cv, fv = setup_vals(ip);
dh = setup_dh(grid, ip);
cst = setup_bcs(dh);

# Extract dofs on the "top" boundary
top = getfacetset(grid, "top")
dofsDtN = dofs_on_dtn(dh, :u, top)

# Allocate A₀, A₁, A₂
A₀ = allocate_stiff_matrix(dh, cst, dofsDtN)
A₁ = allocate_stiff_matrix(dh, cst, dofsDtN)
A₂ = allocate_stiff_matrix(dh, cst, dofsDtN)

# Assemble A₀, A₁, A₂
A₀ = assemble_A0(cv, dh, A₀, medium, k)
A₁ = assemble_A1(cv, dh, A₁)
A₂ = assemble_A2(cv, dh, A₂)

# Impose the boundary conditions
apply!(A₀, cst)
apply!(A₁, cst)
apply!(A₂, cst)

d = size(A₀, 1)

# Construct the nonlinear eigenvalue problem
L = Nnep(A₀, A₁, A₂, fv, dh, cst, top, dofsDtN, N, k);

## Define the contour 
# For the example 1
elp = Cim.circle([0.36, 0.0], 0.05)

# For the example 2
# exceptional value ≈ 0.3803
# elp = Cim.circle([0.38, 0.0], 0.02)

# exceptional value ≈ 0.03
# elp = Cim.circle([1.03, 0.0], 0.05)

# Solve the nonlinear eigenvalue problem by the Contour integral method
λ = new_cim(elp, L, d, 5; n=50, tol=0.5)