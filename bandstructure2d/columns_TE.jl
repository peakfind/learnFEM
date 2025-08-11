# A square lattice of dielectric columns TE modes
# Reference: Photonic Crystals Molding the Flow of Light 2nd, Chapter 5, Figure 2 

using Gmsh
using Ferrite, FerriteGmsh
using CairoMakie
using LinearAlgebra # for norm
using SparseArrays
using Arpack

include("functions.jl")

function ϵ(x)
    r = sqrt(x[1]^2 + (x[2] - 0.5)^2)
    radius = 2π * 0.2

    if r <= radius
        return 8.9
    else
        return 1.0
    end
end

sl = SquareLattice(2π, 2π)
ibz = IrreducibleBrillouin(sl)
dibz, para = get_discrete_IrreducibleBrillouin(ibz, 10)

# Generate the mesh for the periodic cell
grid = setup_grid_squareLattice(lc=0.1, period=2π)

# Add new facetsets for the interior of boundaries
addfacetset!(grid, "itop", is_interior_top)
addfacetset!(grid, "ibottom", is_interior_bottom)
addfacetset!(grid, "ileft", is_interior_left)
addfacetset!(grid, "iright", is_interior_right)

# Get the FacetIndexs corresponding to the new facetsets
itop = getfacetset(grid, "itop")
ibottom = getfacetset(grid, "ibottom")
ileft = getfacetset(grid, "ileft")
iright = getfacetset(grid, "iright")

# Define the interpolation
ip = Lagrange{RefTriangle, 1}()

# Define the CellValues
qr = QuadratureRule{RefTriangle}(2)
cv = CellValues(qr, ip)

# Define the DofHandler
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

# !!! warning 
#     This part is just an Ad hoc implementation. And it only 
#     works for P1 element. 
#     TODO: We need to find a more general way.

# Get the FacetIndexs of the boundaries and the corresponding DoFs
top = getfacetset(grid, "top")
right = getfacetset(grid, "right")
left = getfacetset(grid, "left")
bottom = getfacetset(grid, "bottom")

t = dofs_on_facetset(dh, :u, top) 
r = dofs_on_facetset(dh, :u, right) 
l = dofs_on_facetset(dh, :u, left) 
b = dofs_on_facetset(dh, :u, bottom) 

# Get the DoFs of the four corner points
tr = only(intersect(t, r))
tl = only(intersect(t, l))
br = only(intersect(b, r))
bl = only(intersect(b, l))

# Set up the bi-periodic boundary condition
# 1. Construct the ConstraintHandler
cst = ConstraintHandler(dh)

# 2. Collect the constraints between "iright" and "ileft"
pfacets_x = collect_periodic_facets(dh.grid, "iright", "ileft", x -> x + Ferrite.Vec{2}((2π, 0.0)))
pbc_x = PeriodicDirichlet(:u, pfacets_x)
add!(cst, pbc_x)

# 3. Collect the constraints between "itop" and "ibottom"
pfacets_y = collect_periodic_facets(dh.grid, "itop", "ibottom", x -> x + Ferrite.Vec{2}((0.0, 2π)))
pbc_y = PeriodicDirichlet(:u, pfacets_y)
add!(cst, pbc_y)

# 4. Handle the constraints between the corner points
lc1 = Ferrite.AffineConstraint(tr, [bl => 1.0], 0.0)
lc2 = Ferrite.AffineConstraint(br, [bl => 1.0], 0.0)
lc3 = Ferrite.AffineConstraint(tl, [bl => 1.0], 0.0)
add!(cst, lc1)
add!(cst, lc2)
add!(cst, lc3)

# 5. Close the ConstraintHandler
close!(cst)

# Allocate the matrices 
A = allocate_matries(dh, cst)
B = allocate_matries(dh, cst)

μ = calc_diagram_TE(cv, dh, cst, A, B, ϵ, dibz, nevs=4)
μ = abs.(μ)
μ = sqrt.(μ)
plot_diagram(para, μ; title = "TE mode")