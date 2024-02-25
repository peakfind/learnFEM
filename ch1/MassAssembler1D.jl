"""
    MassAssembler1D(x)

    assemble global mass matrix
"""
function MassAssembler1D(x)
    n = length(x) - 1
    M = zeros(n+1, n+1)

    for i = 1:n # loop over subintervals
        h = x[i+1] - x[i] # interval length
        M[i, i] = M[i, i] + h/3
        M[i, i+1] = M[i, i+1] + h/6 
        M[i+1, i] = M[i+1, i] + h/6
        M[i+1, i+1] = M[i+1, i+1] + h/3 
    end

    return M
end