using Plots
include("L2Projector1D.jl")

f(y) = y*sin(y)
x = range(0, 1, length=100)
y = @. f(x)

nodes = range(0, 1, 6)
l2Proj = L2Projector1D(f)

plot(x, y)
plot!(nodes, l2Proj)

