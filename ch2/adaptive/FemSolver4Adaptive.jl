using LinearAlgebra

function StiffnessAssembler(x)
    n = length(x) - 1
    A = zeros(n+1, n+1)

    for i = 1:n
        h = x[i+1] - x[i]

        A[i, i] = A[i, i] + 1/h
        A[i, i+1] = A[i, i+1] - 1/h 
        A[i+1, i] = A[i+1, i] - 1/h 
        A[i+1, i+1] = A[i+1, i+1] + 1/h
    end

    # apply the zero Dirichlet boundary condition
    A[1,1] = 1.0
    A[n+1,n+1] = 1.0
    A[1,2:end] .= 0.0
    A[n+1,1:n] .= 0.0
    A[2:n,1] .= 0.0
    A[2:n,n+1] .= 0.0

    return A
end

function SourceAssembler(x, f::Function)
    n = length(x) - 1
    b = zeros(n+1)

    for i = 1:n
        h = x[i+1] - x[i]
        b[i] = b[i] + f(x[i])*h/2
        b[i+1] = b[i+1] + f(x[i+1])*h/2
    end

    # apply the zero Dirichlet boundary condition
    b[1] = 0.0
    b[end] = 0.0

    return b
end

function FemSolver(x, f::Function)
    A = StiffnessAssembler(x)
    b = SourceAssembler(x, f)
   
    u = A\b
    return u
end