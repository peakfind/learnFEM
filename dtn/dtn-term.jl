# dtn-term.jl
# Construct a sparse vector Θⁿ
 
using Ferrite, SparseArrays

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
qr = QuadratureRule{RefTriangle}(4)
qr_facet = FacetQuadratureRule{RefTriangle}(2)
cv = CellValues(qr, ip)
fv = FacetValues(qr_facet, ip)

# Define a DofHandler
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

# Extract the dofs which are associated to the DtN map
dofsDtN = dofs_on_dtn(dh, :u, getfacetset(grid, "top"))
Θ = sparsevec(dofsDtN, zeros(ComplexF64, length(dofsDtN)))

# Compute Θⁿ for n = 1
# Allocate the local vector θ
n = getnbasefunctions(fv)
θ = zeros(ComplexF64, n)
for facet in FacetIterator(dh, getfacetset(grid, "top"))
    # Update the fv to the correct facet
    reinit!(fv, facet)
    # Reset the local vector θ
    fill!(θ, 0.0 + 0.0im)
    # 
    coords = getcoordinates(facet)

    # Loop over quadrature points
    for qp in 1:getnquadpoints(fv)
        ds = getdetJdV(fv, qp)
        coords_qp = spatial_coordinate(fv, qp, coords) # Coordinate of the quadrature point
        mode = exp(im*coords_qp[1])
        for i in 1:getnbasefunctions(fv)
            ϕ = shape_value(fv, qp, i)
            θ[i] += ϕ * mode * ds
        end
    end

    assemble!(Θ, celldofs(facet), θ)
    @show celldofs(facet)
    @show θ
end

@show Θ