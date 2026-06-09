using CairoMakie
using NumericalSteepestDescent
using Polynomials

function InterpolatedPhase(a,b,g; n_interp)
    # x,_ = gausschebyshevt(n_interp)
    x = [cos((k+1/2)*π/n_interp) for k in 0:n_interp-1]
    xi = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    P = fit(Polynomial, xi, g.(xi))
    @show degree(P)
    return PolynomialPhase(P.coeffs)
end

g(z) = 1/(z^2 + 1)
gi = InterpolatedPhase(-1,1,g; n_interp = 12)

fl = NumericalSteepestDescent.plot_landscape(gi)
limits!(-2,2,-2,2)
fl

NumericalSteepestDescent.nsd([-1, 1], x -> 1, gi, ω; N = 20, 
                        infquadrule = :lag, plot_sd = true)


ε = 0.01
n = 4
g(z) = -ε*z^(n+2)/(n+2) + z^(n+1)/(n+1) + ε*z^2/2 -z 
gi = InterpolatedPhase(-1,1,g; n_interp = n+3)
println("degree: ", degree(gi.p))
println("largest stationary point: ", 1/ε)
println("rstar: ", gi.rstar_valley)

""" OBSERVATIONS 
- rstar_valley = O(1/ε)
- perhaps we want to keep distance of stat points to origin reasonably similar?

"""


fl = NumericalSteepestDescent.plot_landscape(gi; ran = 2/ε, color_lim = 1e3)
# θ = range(0,2π, length = 100)
# rstar = gi.rstar_valley
# lines!(reim.(rstar*cis.(θ)), color = :red)
# limits!(-2,2,-2,2)
fl

ω = 10
i1, f1 = NumericalSteepestDescent.nsd([-1, 101], x -> 1, gi, ω; N = 20, infquadrule = :lag, plot_sd = true)
f1[1]
limits!(-5,12,-5,5)
f1[1]