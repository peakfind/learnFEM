# Test the code snippets for DtN-FEM in dtn.jl
# An flat surface illuminated by an plane wave with 
# Dirichlet boundary condition.

using Ferrite
using SparseArrays
using Gmsh
using FerriteGmsh

include("dtn.jl")

## Parameters
# Incident field
k = 1.0
θ = π/3
inc = Incident(k, θ)

# Number of the truncated terms
N = 10
height = 2.0
d = 2π

# Generate a uniform mesh
grid = periodic_cell(0.1; height=height)
# grid = generate_grid(Triangle, (40, 40), Vec{2}((0.0, 0.0)), Vec{2}((d, 0.0)), Vec{2}((d, height)), Vec{2}((0.0, height)))

## Set up fevalues(CellValues and FacetValues), DofHandler, and ConstraintHandler
ip = Lagrange{RefTriangle, 1}() # Interpolation and quadrature rules
cv, fv = setup_vals(ip) 
dh = setup_dofs(grid, ip)
cst = setup_bcs(dh)

## Allocate the stiffness matrix A, the TBC matrix F and the load vector f 
# Extract dofs on the "top" boundary
dofsDtN = dofs_on_dtn(dh, :u, getfacetset(grid, "top"))

# Allocate the stiffness matrix and the load vector
A = allocate_stiff_matrix(dh, cst, dofsDtN)
F = allocate_stiff_matrix(dh, cst, dofsDtN)
f = zeros(ComplexF64, ndofs(dh))

# Assemble the matrix A
A = assemble_A(cv, dh, A, inc)

# Assemble the load vector f
f = assemble_load(fv, dh, getfacetset(grid, "top"), f, inc, height)

# Assemble the TBC matrix
F = assemble_tbc(fv, dh, inc, getfacetset(grid, "top"), F, N, dofsDtN)

A = sub_preserve_structure(A, F)

# Get the global stiffness matrix by A - F
# A .-= F

# Impose the boundary condition
apply!(A, f, cst)

# Solve the linear system
u = A\f 
apply!(u, cst)

# Write the solution to vtk file
VTKGridFile("real_u", grid) do vtk
   write_solution(vtk, dh, real.(u)) 
end;

VTKGridFile("imag_u", grid) do vtk
    write_solution(vtk, dh, imag.(u)) 
end;