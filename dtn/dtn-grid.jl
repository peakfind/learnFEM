#=
# dtn-grid.jl 
#
# Learn the basic notions in Ferrite.jl such as gird, dofs. ...
=#

using Ferrite, FerriteGmsh, Gmsh, SparseArrays

function setup_grid(lc=0.5)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(2, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(2, 2, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, 2, 0, lc)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom") 
    gmsh.model.addPhysicalGroup(1, [l2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], [1, 0, 0, 2, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.generate(2)
    grid = mktempdir() do dir
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path) 
    end
    gmsh.finalize()
    return grid
end

function setup_vals(interpolation)
    qr = QuadratureRule{RefTriangle}(2)
    qr_facet = FacetQuadratureRule{RefTriangle}(2)
    cellvalues = CellValues(qr, interpolation)
    facetvalues = FacetValues(qr_facet, interpolation)
    return cellvalues, facetvalues
end

function setup_dofs(grid::Grid, interpolation)
    dh = DofHandler(grid)
    add!(dh, :u, interpolation)
    close!(dh)
    return dh
end

function setup_bcs(dofhandler::DofHandler) 
    ch = ConstraintHandler(dofhandler)
    periodic_faces = collect_periodic_facets(dofhandler.grid, "right", "left", x -> x + Vec{2}((2, 0.0)))
    pbc = PeriodicDirichlet(:u, periodic_faces)
    add!(ch, pbc)
    dbc = Dirichlet(:u, getfacetset(dofhandler.grid, "bottom"), x -> 0)
    add!(ch, dbc)
    close!(ch)
    return ch
end

function dofs_on_dtn(dofhandler::DofHandler, field::Symbol, facetset)
    dtn_ch = ConstraintHandler(dofhandler)
    dbc = Dirichlet(field, facetset, x -> 0)
    add!(dtn_ch, dbc)
    close!(dtn_ch)
    return dtn_ch.prescribed_dofs
end

function allocate_stiff_matrix(dofhandler::DofHandler, csthandler::ConstraintHandler, dofs)
    sp = init_sparsity_pattern(dofhandler)
    add_cell_entries!(sp, dofhandler)
    # Use add_entry! for DtN term
    for i in dofs, j in dofs
        if abs(i - j) > 1
            Ferrite.add_entry!(sp, i, j)
        end
    end
    add_constraint_entries!(sp, csthandler)
    K = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return K
end

function main()
    grid = setup_grid()
    ip = Lagrange{RefTriangle, 1}()
    dh = setup_dofs(grid, ip)
    ch = setup_bcs(dh)
    dofs_dtn = dofs_on_dtn(dh, :u, getfacetset(grid, "top"))
    K = allocate_stiff_matrix(dh, ch, dofs_dtn)
    return K
end

main()