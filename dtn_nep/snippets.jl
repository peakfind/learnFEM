#----------------------------------------------------------
# Finite element discretization
#----------------------------------------------------------

"""
    assemble_A0(cv::CellValues, dh::DofHandler, A₀::SparseMatrixCSC, medium::Function, k)

Assemble the zero order term ``\\mathbf{A}_{0}`` in the Nonlinear eigenvalue problem.

``(\\mathbf{A}_{0} + \\alpha \\mathbf{A}_{1} + \\alpha^2 \\mathbf{A}_{2} + \\mathbf{F}(\\alpha)) \\mathbf{u} = 0.``
"""
function assemble_A0(cv::CellValues, dh::DofHandler, A₀::SparseMatrixCSC, medium::Function, k)
    # Allocate the local matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)
    
    # Create an assembler A₀ 
    assembler = start_assemble(A₀)
    
    # Loop over all cells
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)

        # Reset local stiffness matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)
        coords = getcoordinates(cell)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            n = medium(coords_qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                ∇v = shape_gradient(cv, qp, i)

                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble local stiffness matrix
                    Ae[i, j] += (∇u[1] * ∇v[1] + ∇u[2] * ∇v[1] - (k^2) * n * u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₀ 
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₀
end

"""
    assemble_A1(cv::CellValues, dh::DofHandler, A₁::SparseMatrixCSC)

Assemble the first order term ``\\mathbf{A}_{1}`` in the Nonlinear eigenvalue problem.

`` (\\mathbf{A}_{0} + \\alpha \\mathbf{A}_{1} + \\alpha^2 \\mathbf{A}_{2} + \\mathbf{F}(\\alpha)) \\mathbf{u} = 0. ``
"""
function assemble_A1(cv::CellValues, dh::DofHandler, A₁::SparseMatrixCSC)
    # Allocate the local stiffness matrix 
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)

    # Create an assembler 
    assembler = start_assemble(A₁)
    
    # Loop over all cells
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        
        # Reset local stiffness matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)
        
        # Loop over quadrature points 
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs 
                v = shape_value(cv, qp, i)

                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble local stiffness matrix 
                    Ae[i, j] += (-2im * ∇u[1] * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₁
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₁
end

"""
    assemble_A2(cv::CellValues, dh::DofHandler, A₂::SparseMatrixCSC)

Assemble the second order term in the Nonlinear eigenvalue problem.

`` (\\mathbf{A}_{0} + \\alpha \\mathbf{A}_{1} + \\alpha^2 \\mathbf{A}_{2} + \\mathbf{F}(\\alpha)) \\mathbf{u} = 0. ``
"""
function assemble_A2(cv::CellValues, dh::DofHandler, A₂::SparseMatrixCSC)
    # Allocate the local matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)

    # Create an assembler A₂
    assembler = start_assemble(A₂)

    # Loop over all cells 
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell 
        reinit!(cv, cell)
        
        # Reset local stiffness matrix to 0.0 + 0.0im 
        fill!(Ae, 0.0 + 0.0im)
        
        # Loop over quadrature points 
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                
                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    
                    # Assemble local stiffness matrix 
                    Ae[i, j] += (u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₂ 
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₂
end

#----------------------------------------------------------
# Nonlinear eigenvalue problem
#----------------------------------------------------------
"""
    Nep

Nep is used to store the matrices assocaited to the quadratic part.

    (A₀ + α A₁ + α² A₂ + F(α))u = 0 
"""
struct Nep 
    A₀::SparseMatrixCSC
    A₁::SparseMatrixCSC
    A₂::SparseMatrixCSC
end

# TODO: construct nonlinear eigenvalue problem u::Nep by u = Nep(A0, A1, A2)
# TODO: function (L::Nep)(\alpha::Complex)

"""
    (nep::Nep)(α)

TBW
"""
function (nep::Nep)(α)
    # TODO:Generate TBC matrix F(α)

    N = nep.A₀ + α * A₁ + (α^2) * A₂
    return N
end