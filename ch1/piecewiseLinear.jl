using Plots

x = range(0, 1, length=100)
nodes = range(0, 1, length=4)

## f = z^2 + 1
# f(z) = z*z + 1
# y = @. f(x)
# y_i = @. f(nodes)

## f = cos(pi*z)
f(z) = cos(pi*z)
y = @. f(x)
y_i = @. f(nodes)

plot(x, y)
plot!(nodes, y_i)