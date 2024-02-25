using Plots

x = range(0, 1, length=100)

## y = x^2
# y = @. x*x
# y_i = x

## y = 3 sin(2\pi x)
y = @. 3*sin(2*pi*x)
y_i = 0*x

plot(x, [y, y_i])