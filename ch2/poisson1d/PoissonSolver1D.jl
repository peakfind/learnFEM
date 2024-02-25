using LinearAlgebra
using Plots

include("StiffnessAssembler1D.jl")
include("SourceAssembler1D.jl")

function PoissonSolver1D()
    x = range(2, 8, length=26)
    Conductivity(x) = 0.1*(5 - 0.6*x)
    Source(x) = 0.03*(x-6)^4
    kappa = [1.e+6, 0]
    g = [-1, 0]

    A = StiffnessAssembler1D(x, Conductivity, kappa)
    b = SourceAssembler1D(x, Source, kappa, g)
    u = A\b

    plot(x, u)
end