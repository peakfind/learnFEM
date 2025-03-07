using Gmsh
using Ferrite, FerriteGmsh
using SparseArrays

include("periodic_surface.jl")

# Parameters
k = 2.0
θ = 0.19π

# Setup the grid
grid = setup_grid(0.1)

# Interpolation
ip = Lagrange{RefTriangle, 1}()

# Set the FE values
cv = setup_fevs(ip)

# Set the DofHandlers
dh = setup_dofs(grid, ip)

# Set boundary conditions
cst = setup_bdcs(dh, k, θ)

# Allocate the stiffness matrix and right hand side term
A = allocate_matrices(dh, cst)
f = zeros(ComplexF64, ndofs(dh))

A = assemble_global!(cv, dh, A, k, θ)
apply!(A, f, cst)
u = A\f 
apply!(u, cst)

VTKGridFile("real_u", grid) do vtk
    write_solution(vtk, dh, real(u))
end;

VTKGridFile("imag_u", grid) do vtk
    write_solution(vtk, dh, imag(u))
end;

VTKGridFile("modulus_u", grid) do vtk
    write_solution(vtk, dh, abs.(u))
end;