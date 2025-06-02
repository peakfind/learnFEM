using Ferrite
# reuse some functions
using Nmrc: periodic_cell, setup_vals, setup_dh, setup_bcs, dofs_on_dtn, allocate_stiff_matrix,
            assemble_A0, assemble_A1, assemble_A2, compute_coef!
using Cim
using LinearAlgebra
using SparseArrays
using Random

include("snippets.jl")

# Wavenumber
k = 4.1

# The nb of the truncated terms
N = 10

# Refractive index
function medium(x)
    if x[2] < 1.0
        return 3.9
    else
        return 1.0
    end
end

# Generate the mesh in a periodic cell
grid = periodic_cell(lc=0.03, period=2π, height=2.0)

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

# Define the contour 
elp = Cim.circle([0.38, 0.0], 0.05)

# Solve the nonlinear eigenvalue problem by the Contour
# integral method
λ = new_cim(elp, L, d, 5; n=200, tol=0.9)
@show λ