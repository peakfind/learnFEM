using Plots
include("FemSolver4Adaptive.jl")

# refinement parameter
alpha = 0.9

# the source term
c = 100
Source(z) = exp(-c*(z - 0.5)^2)

# initial mesh
x = range(0, 1, 5) 

MaxIter = 25

for iter = 1:MaxIter
    global x
    n = length(x) - 1
    eta = zeros(n) # allocate element residual
    for i = 1:n
        h = x[i+1] - x[i]
        t = (Source(x[i])^2 + Source(x[i+1])^2)*h/2
        eta[i] = h*sqrt(t)
    end
    for i in eachindex(eta)
        if eta[i] > alpha*maximum(eta)
            # insert new node
            x = vcat(x,[(x[i+1] + x[i])/2])
        end
    end
    # get the new mesh
    x = sort(x) 
end

uAdaptive = FemSolver(x, Source)
scatter!(x, uAdaptive, mc=:red)

