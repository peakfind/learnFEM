# dtn-grid.jl
# Goal: Extract dofs on a specific facetset
# References: [boundary dof #555](https://github.com/Ferrite-FEM/Ferrite.jl/discussions/555)

using Ferrite

function dofs_on_dtn(dofhandler::DofHandler, field::Symbol, facetset)
    dtn_ch = ConstraintHandler(dofhandler)
    dbc = Dirichlet(field, facetset, x -> 0)
    add!(dtn_ch, dbc)
    close!(dtn_ch)
    return dtn_ch.prescribed_dofs
end

# Generate a uniform mesh
grid = generate_grid(Triangle, (10, 10))

# Combine the interpolation and the quadrature rule to CellValues
ip = Lagrange{RefTriangle, 1}()
qr = QuadratureRule{RefTriangle}(2)
cv = CellValues(qr, ip)

# Define a DofHandler
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

dofsDtN = dofs_on_dtn(dh, :u, getfacetset(grid, "top"))
for facet in FacetIterator(dh, getfacetset(grid, "top"))
     @show getcoordinates(facet)
end