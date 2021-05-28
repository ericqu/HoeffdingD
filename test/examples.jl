# Readme example 

x = -2:0.1:2
linear_f(x) = 2x ; quad_f(x) = x^2
y_linear = linear_f.(x)
y_quad = quad_f.(x)

using Plots
scatter(x, y_linear, label="linear")
scatter!(x, y_quad, label="quadratic", legend=:bottomright)
savefig("docs/linquad.png")

using StatsBase

@show(StatsBase.cor(x, y_linear))
@show(StatsBase.cor(x, y_quad))
@show(StatsBase.corspearman(x, y_linear))
@show(StatsBase.corspearman(x, y_quad))


@show(HoeffdingD.hoeffdingd(x, y_linear))
@show(HoeffdingD.hoeffdingd(x, y_quad))
@show(HoeffdingD.hoeffdingd(x, y_linear, 0.05))
@show(HoeffdingD.hoeffdingd(x, y_quad, 0.05))

