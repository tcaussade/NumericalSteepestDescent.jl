using Polynomials
using FastGaussQuadrature
using CairoMakie
using QuadGK
using NumericalSteepestDescent


function InterpolatedPhase(a,b,g; n_interp)
    x,_ = gausschebyshevt(n_interp)
    xi = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    P = fit(Polynomial, xi, g.(xi))
    @show degree(P)
    return PolynomialPhaseFunction(P.coeffs)
end

# x,_ = gausschebyshevt(20)
# g(z) = cos(z)
# P1 = fit(ChebyshevT, x, g.(x))
# P2 = fit(Polynomial, x, g.(x))

# idx = P2.coeffs .> 1e-6
# P2t = P2[idx]

function plot_nsd(a,b,g :: NumericalSteepestDescent.AbstractPhaseFunction, ω)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Re(z)", ylabel = "Im(z)", title = "NSD for polynomial interpolation", aspect = 1)
    quasiSDdeformation!(fig, ax, [a,b], g, ω; 
                        umax = 10, resolution = 200, color_lim = 100)
    # limits!(-5,5,-5,5)
    # scatter!(ax, reim.(-5π:π:5π), color = :red, markersize = 12, marker = :star5)
    scatter!(ax, reim.([-im,im]), color = :red, markersize = 12, marker = :star5)
    return fig
end

function plot_landscape(g, η)
    nx = 100
    color_lim = 100
    x = range(-10,10, length = nx)
    y = range(-10,10, length = nx)
    X = [x for x in x for _ in y]
    Y = [y for _ in x for y in y]
    Z = [g.p(x+im*y) for x in x for y in y]

    rlevels = [unique(real.(g.p.(g.ξ))); real(g.p(η))]


    # color_lim = 300 # maximum(imag(Z)) / 2
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Re(z)", ylabel = "Im(z)", title = "SD landscape", aspect = 1)
    levelset = contourf!(ax,X,Y,-imag.(Z); levels = range(-color_lim, color_lim, 20), 
                         colormap = :balance, extendlow = :auto, extendhigh = :auto)
    contour!(ax,X,Y,real.(Z); levels = rlevels, color = :black, linewidth = 1, linestyle = :dash)
    scatter!(ax, reim.(g.ξ), color = :red, markersize = 12, marker = :star5) 
    scatter!(ax, reim.([η]), color = :blue, markersize = 12, marker = :circle)
    Colorbar(fig[1,2], levelset)
    return fig
end


# Basic test in [a,b] for g(z) = cos(z)
a,b = -1.,+1.0
ω = 30
g(z) = 1/(1+z^2)
# g(z) = cos(z)
g1 = InterpolatedPhase(a,b,g; n_interp = 52)

fl = plot_landscape(g1,1)
limits!(-2,2,-2,2)
fl

igk = quadgk(x -> cis(ω * g(x)), a, b, rtol = 1e-14)[1]
i1 = NumericalSteepestDescent.integrate([a, b], x -> 1, g1, ω; N = 20, infquadrule = :lag,
    δODE = 0.01)
@show abs(i1 - igk)/abs(igk)

f1 = plot_nsd(a,b,g1,ω)
limits!(-2,2,-2,2)
f1

"""
Observations:
- SD contour tracing may fail when αJ is small, because rstar_valley becomes very large. 
This causes newton iteration to fail before reaching the valley.
- For large polynomial degree, the computation of roots to solve the ball radius becomes unstable, 
and the "real root" that we are looking for may have large imaginary part. 
"""

# Test against DLMF 10.9.2
using SpecialFunctions

ν = 0
z = 15
ref = besselj(ν, z)

g2 = InterpolatedPhase(0,π,g; n_interp = 14)
i2= im^(-ν)/π * NumericalSteepestDescent.integrate([0,π], x -> cos(ν*x), g2, z)
@show abs(i2 - ref)/abs(ref)s

f2 = plot_nsd(0,π,g2,z)
display(f2)

""" TESTS TO DO:
- plot laguerre nodes and see where do they go along SD contour
- Test error vs. n_interp (keep Nquad large)
- Test error vs. Nquad (keep n_interp fixed)
- Test error vs. ω (keep n_interp and Nquad fixed)
Repeat for g(z)=cos(z) and g(z)=1/(1+z^2)
"""
