using LinearAlgebra
include("MassAssembler1D.jl")
include("LoadAssembler1D.jl")

function L2Projector1D(f::Function)
    n = 5
    x = LinRange(0, 1, n+1)
    M = MassAssembler1D(x)
    b = LoadAssembler1D(x, f)
    Pf = M\b

    return Pf
end