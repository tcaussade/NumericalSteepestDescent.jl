
""" TESTS TO DO:
- plot laguerre nodes and see where do they go along SD contour
- Test error vs. n_interp (keep Nquad large)
- Test error vs. Nquad (keep n_interp fixed)
- Test error vs. ω (keep n_interp and Nquad fixed)
Repeat for g(z)=cos(z) and g(z)=1/(1+z^2)
"""

using Polynomials
using FastGaussQuadrature
using CairoMakie
using QuadGK
using NumericalSteepestDescent

function InterpolatedPhase(a,b,g; n_interp)
    x,_ = gausslegendre(n_interp)
    xi = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    P = fit(Polynomial, xi, g.(xi))
    # @show degree(P)
    return PolynomialPhaseFunction(P.coeffs)
end

function experiment_interp(x,g, nquad, ninterp, ω; sdquad)
    # endpoints are [a,b]
    # it is assumed that a < x1 < x2 < ... < xn < b
    # x = [a,x1,x2,...,xn,b]
    gval = zero(ComplexF64)
    for i in eachindex(x)[1:end-1]
        x1,x2 = x[i], x[i+1]
        @show "doing" x1,x2
        gint = InterpolatedPhase(x1,x2,g; n_interp = ninterp)
        gval += NumericalSteepestDescent.integrate([x1, x2], x -> 1, gint, ω; N = nquad, infquadrule = sdquad)
    end
    return gval
end

# Test for functions
a,b = -1,1
g(z) = 1/(1+z^2)
# g(z) = cos(z)
# g(z) = sqrt(z+2)

# Test interpolation errors
nquad = 20
ωvals = [5.0, 10.0, 20.0, 40.]
cols = [:blue, :orange, :green, :red]
ninterp_vals = 5:2:20

fig = Figure()
ax1 = Axis(fig[1, 1], xlabel = "n_interp", ylabel = "absolute error", title = "Error vs. n_interp",
    yscale = log10)

x = [a,b]
xvals = ninterp_vals * (length(x)-1)
for (i,ω) in enumerate(ωvals)
    vals_lag = [experiment_interp(x, g, nquad, ninterp, ω; sdquad = :lag) for ninterp in ninterp_vals]
    vals_tleg = [experiment_interp(x, g, nquad, ninterp, ω; sdquad = :tleg) for ninterp in ninterp_vals]
    ref = quadgk(x -> cis(ω * g(x)), a, b, rtol = 1e-14)[1]
    errors_lag = abs.(vals_lag .- ref) ./ abs.(ref)
    errors_tleg = abs.(vals_tleg .- ref) ./ abs.(ref)

    scatterlines!(ax1, xvals, errors_lag, label = "ω = $ω", color = cols[i], linestyle = :solid)
    scatterlines!(ax1, xvals, errors_tleg, color = cols[i], linestyle = :dash)
end
axislegend(ax1) 
limits!(2, 21, 1e-16, 1e-0)
fig

# Test subdivision of intervals e.g. with equispacing
ndiv = [1, 2, 4]
ninterp = 5

fig = Figure()
ax1 = Axis(fig[1, 1], xlabel = "sum of ninterps", ylabel = "absolute error", title = "Error vs. n_interp",
    yscale = log10)
for (i,ω) in enumerate(ωvals)
    ref = quadgk(x -> cis(ω * g(x)), a, b, rtol = 1e-14)[1]
    ndiv_error = Float64[]
    for n in ndiv
        @show ω,n
        x = range(a, b, length = n+1)
        gval = experiment_interp(x, g, nquad, ninterp, ω; sdquad = :lag)
        push!(ndiv_error, abs(gval .- ref)./abs(ref))
    end
    scatterlines!(ax1, ndiv * ninterp, ndiv_error, label = "ω = $ω", color = cols[i], linestyle = :solid)
end
axislegend(ax1, position = :lb) 
limits!(2, 41, 1e-16, 1e-0)
fig


