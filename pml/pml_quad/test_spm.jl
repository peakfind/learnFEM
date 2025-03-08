using Gmsh
using Ferrite, FerriteGmsh
using SparseArrays
using NonlinearEigenproblems
# using Plots

include("pml_quad.jl")

d = 2π
grid = setup_grid(;d=2π, h=1.5, δ=0.5, lc=0.5)

# Define the interpolation: linear lagrange
ip = Lagrange{RefTriangle, 1}()

# Setup FE values, dofs, and boundary conditions
cv = setup_fevs(ip)
dh = setup_dofs(grid, ip)
ch = setup_bdcs(dh, d)

# Allocate the matrices
f₀ = zeros(ComplexF64, ndofs(dh))
f₁ = zeros(ComplexF64, ndofs(dh))
f₂ = zeros(ComplexF64, ndofs(dh))
A₀ = allocate_matrices(dh, ch)
A₁ = allocate_matrices(dh, ch)
A₂ = allocate_matrices(dh, ch)

A₀ = assemble_A0(cv, dh, A₀)
A₁ = assemble_A1(cv, dh, A₁)
A₂ = assemble_A0(cv, dh, A₂)
# apply!(A₀, f, ch)
# apply!(A₁, f, ch)
# apply!(A₂, f, ch)
apply!(A₀, f₀, ch)
apply!(A₁, f₁, ch)
apply!(A₂, f₂, ch)

# Set NEP 
qep = PEP([A₀, -A₁, A₂])
λ,_ = polyeig(qep)
filter!(e -> is_real(e), λ)
@show λ