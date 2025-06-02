# # Finite element method with transparent boundary condition

# ## Introduction

# ## Numerical examples

# First, we load Nmrc and Ferrite (Ferrite is used to define the interpolation).
using Nmrc, Ferrite

# ### Parameters 

# We assume that the wavenumber is ``1.0`` and the incident angle is $\pi/3$.
# ```math
# \begin{align}
#     x + y \\
#     y + z
# \end{align}
# ```
## incident plane wave
k = 1.0
θ = π/3
inc = Incident(k, θ)

## Number of the truncated terms
N = 10

## Height of the artificial boundary
height = 2.0

# ### Generating a mesh by Gmsh

# We generate a mesh on a periodic cell which is a rectangle $(0.0, 2\pi) \times 
# (0.0, 2.0)$.

## Generate the mesh in a periodic cell
grid = periodic_cell(lc=0.05, period=2π ,height=height)

# ### Setting up finite element space and boundary conditions

# We choose $P1$ element here. Then `CellValues` and `FacetValues` can be constructed 
# by the function [`setup_vals`](@ref).
ip = Lagrange{RefTriangle, 1}() 
cv, fv = setup_vals(ip) 

# `DofHandler` handles numbering and distribution of degrees of freedom (dofs) of 
# our solution.
dh = setup_dh(grid, ip)

# We use [`setup_bcs`](@ref) to construct our `ConstraintHandler` which contains 
# the information of the boundary conditions. In this tutorial, we impose the Dirichelt 
# boundary condition on the lower boundary and the periodic boundary condition on the 
# left and right boundaries.
cst = setup_bcs(dh)

# ### Assembling the stiffness matrix

# We allocate the stiffness matrix $A$, the TBC matrix $F$ and the load vector 
# $f$. Due to the DtN map, we need to extract the dofs associated to the artificial 
# boundary. Then we create the sparse pattern of $A$ and $F$. For more details, 
# please refer to ['allocate_stiff_matrix'](@ref).

## Extract dofs on the "top" boundary
top = getfacetset(grid, "top")
dofsDtN = dofs_on_dtn(dh, :u, top)

## Allocate the stiffness matrix and the load vector
A = allocate_stiff_matrix(dh, cst, dofsDtN)
F = allocate_stiff_matrix(dh, cst, dofsDtN)
f = zeros(ComplexF64, ndofs(dh))

# The matrix $A$ comes from 
# ```math
#   \int dx.
# ```
## Assemble the matrix A
A = assemble_A(cv, dh, A, inc)

## Assemble the load vector f
f = assemble_load(fv, dh, top, f, inc, height)

## Assemble the TBC matrix
F = assemble_tbc(fv, dh, inc, top, F, N, dofsDtN)

# Be careful! Calculations between two sparse matrices may destory the structure 
# of sparse matrix. Then imposing the boundary condition may fail. So we use 
# [`sub_preserve_structure`](@ref) to subtract $F$ from $A$ with keeping the 
# structure of $A$.

## Add the TBC matrix to A = A - F
A = sub_preserve_structure(A, F)

## Impose the boundary condition
apply!(A, f, cst)

## Solve the linear system
u = A\f 
## Apply the boundary condition again
apply!(u, cst)

# ### Export the VTK

# Finally, we write the real part and imaginary part of the solution into vtk 
# files which can be visualized by Paraview.

VTKGridFile("real_u", grid) do vtk
   write_solution(vtk, dh, real.(u)) 
end;

VTKGridFile("imag_u", grid) do vtk
    write_solution(vtk, dh, imag.(u)) 
end;