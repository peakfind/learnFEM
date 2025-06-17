using CairoMakie

## Parameters
# the number of the truncated terms
N = 6

# the wavenumber
k = 4.1

## quadrature nodes on the contour

# the center of the circle (contour)
center = [0.38, 0.0] 
# center = [0.03, 0.0]

# the radius of the circle (contour)
r = 0.02 

num_quadpts = 100
δ = 2π / (num_quadpts - 1)
α = zeros(ComplexF64, num_quadpts)

for i = 0:num_quadpts-1
    α[i + 1] =  complex(center[1] + r * cos(δ * i), center[2] + r * sin(δ * i))
end

fig = Figure()
ax = Axis(fig[1,1], title = L"$\gamma_{n}$ for all the truncated terms", xticks=-100:5:20)

for n in -N:N
    γ = k^2 .- (α .+ n).^2
    scatter!(ax, real.(γ), imag.(γ), color = :red)
end

fig