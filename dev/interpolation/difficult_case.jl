using CairoMakie
using NumericalSteepestDescent

function plot_landscape(g, η; ran = 10, color_lim = 100)
    nx = 200
    x = range(-ran,ran, length = nx)
    y = range(-ran,ran, length = nx)
    X = [x for x in x for _ in y]
    Y = [y for _ in x for y in y]
    Z = [g.p(x+im*y) for x in x for y in y]

    rlevels = [unique(real.(g.p.(g.ξ))); real(g.p(η))]


    # color_lim = 300 # maximum(imag(Z)) / 2
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Re(z)", ylabel = "Im(z)", title = "SD landscape", aspect = 1)
    levelset = contourf!(ax,X,Y,-imag.(Z); levels = range(-color_lim, color_lim, 20), 
                         colormap = :jet, extendlow = :auto, extendhigh = :auto)
    contour!(ax,X,Y,real.(Z); levels = rlevels, color = :black, linewidth = 1, linestyle = :dash)
    scatter!(ax, reim.(g.ξ), color = :red, markersize = 12, marker = :star5) 
    scatter!(ax, reim.([η]), color = :blue, markersize = 12, marker = :circle)
    Colorbar(fig[1,2], levelset)
    return fig
end

function InterpolatedPhase(a,b,g; n_interp)
    # x,_ = gausschebyshevt(n_interp)
    x = [cos((k+1/2)*π/n_interp) for k in 0:n_interp-1]
    xi = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    P = fit(Polynomial, xi, g.(xi))
    @show degree(P)
    return PolynomialPhaseFunction(P.coeffs)
end

g(z) = 1/(z^2 + 1)
gi = InterpolatedPhase(-1,1,g; n_interp = 12)

fl = plot_landscape(gi,1)
limits!(-2,2,-2,2)
fl


ε = 0.001
n = 8
g(z) = -ε*z^(n+2)/(n+2) + z^(n+1)/(n+1) + ε*z^2/2 -z 
gi = InterpolatedPhase(-1,1,g; n_interp = n+3)
println("degree: ", degree(gi.p))
println("largest stationary point: ", 1/ε)
println("rstar: ", gi.rstar_valley)

""" OBSERVATIONS 
- rstar_valley = O(1/ε)
- perhaps we want to keep distance of stat points to origin reasonably similar?

"""


fl = plot_landscape(gi,1; ran = 2/ε, color_lim = 1e3)
θ = range(0,2π, length = 100)
rstar = gi.rstar_valley
lines!(reim.(rstar*cis.(θ)), color = :red)

# limits!(-2,2,-2,2)
fl

ω = 10
i1, f1 = NumericalSteepestDescent.integrate([-1, 11], x -> 1, gi, ω; N = 20, infquadrule = :lag, plot_sd = true)
f1[1]
limits!(-5,12,-5,5)
f1[1]