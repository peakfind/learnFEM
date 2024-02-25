using Plots

function showmesh(P::AbstractArray{Float64,2}, T::AbstractArray{Int64,2})
    nt = size(T, 2)

    # loop in elements
    for k in 1:nt
        plot!(P[1, [T[1,k], T[2,k]]], P[2, [T[1,k], T[2,k]]], legend=:false, showaxis=:false, lc=:black, show=:true, grid=:false)
        plot!(P[1, [T[2,k], T[3,k]]], P[2, [T[2,k], T[3,k]]], legend=:false, showaxis=:false, lc=:black, show=:true, grid=:false)
        plot!(P[1, [T[3,k], T[1,k]]], P[2, [T[3,k], T[1,k]]], legend=:false, showaxis=:false, lc=:black, show=:true, grid=:false)
    end
end