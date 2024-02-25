using LinearAlgebra
using Plots

"""
    Here are the three functions used for L2 projection. We 
    introduce some notations that we will use in the codes:
    
    scalars
    -------------------------------------------------------
    np: the number of nodes
    nt: the number of elements

    1d arrary
    -------------------------------------------------------
    ltg: the local-to-global mapping between local numbers 
         and global numbers

    2d arrary
    -------------------------------------------------------
    P: the point matrix P is of size 2 x np and contains the 
       coordinates of all nodes
    T: the connectivity matrix T is of size 3 x nt and contains 
       the index of the three nodes in each triangle
    
"""

function MassAssembler2D(P::AbstractArray{Float64,2}, T::AbstractArray{Int64,2})
    np = size(P, 2)
    nt = size(T, 2)

    # allocate mass matrix TODO:use sparse matrix in julia?
    M = zeros(np, np)

    # loop over elements
    for k in 1:nt
        # set up the local-to-global mapping
        ltg = T[:,k]

        # get the coordinates of all nodes in triangle k
        x1 = P[1, ltg]
        x2 = P[2, ltg]

        # compute the area of the triangle k
        area = abs((x1[2] - x1[1])*(x2[3] - x2[1]) - (x2[2] - x2[1])*(x1[3] - x1[1]))/2

        # add the local element mass matrix to the global mass matrix
        M[ltg[1], ltg[1]] = M[ltg[1], ltg[1]] + area/6
        M[ltg[1], ltg[2]] = M[ltg[1], ltg[2]] + area/12
        M[ltg[1], ltg[3]] = M[ltg[1], ltg[3]] + area/12

        M[ltg[2], ltg[1]] = M[ltg[2], ltg[1]] + area/12
        M[ltg[2], ltg[2]] = M[ltg[2], ltg[2]] + area/6
        M[ltg[2], ltg[3]] = M[ltg[2], ltg[3]] + area/12

        M[ltg[3], ltg[1]] = M[ltg[3], ltg[1]] + area/12
        M[ltg[3], ltg[2]] = M[ltg[3], ltg[2]] + area/12
        M[ltg[3], ltg[3]] = M[ltg[3], ltg[3]] + area/6
    end

    return M
end 

function LoadAssembler2D(P::AbstractArray{Float64,2}, T::AbstractArray{Int64,2}, f::Function)
    np = size(P, 2)
    nt = size(T, 2)

    # allocate the load vector
    b = zeros(np)

    # loop over elements
    for k in 1:nt
        # set up the local-to-global mapping
        ltg = T[:,k]

        # get the coordinates of all nodes in triangle k
        x1 = P[1, ltg]
        x2 = P[2, ltg]

        # compute the area of the triangle k
        area = abs((x1[2] - x1[1])*(x2[3] - x2[1]) - (x2[2] - x2[1])*(x1[3] - x1[1]))/2

        # add the local element load vector to the global load vector
        b[ltg[1]] = b[ltg[1]] + f(x1[1], x2[1])*area/3
        b[ltg[2]] = b[ltg[2]] + f(x1[2], x2[2])*area/3
        b[ltg[3]] = b[ltg[3]] + f(x1[3], x2[3])*area/3
    end

    return b
end

function L2Projector2D(P::AbstractArray{Float64,2}, T::AbstractArray{Int64,2}, f::Function)
    M = MassAssembler2D(P, T)
    b = LoadAssembler2D(P, T, f)

    L2P = M\b
end