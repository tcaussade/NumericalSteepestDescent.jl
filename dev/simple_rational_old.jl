using NumericalSteepestDescent
using QuadGK
using FastGaussQuadrature
using CairoMakie

ω = 2.

# trying g(z) = z + 1/z
acoefs = [0,1]
ps     = [0.0]
pcoefs = [[1]] #

RatPhase = RationalPhaseFunction(acoefs, ps, pcoefs)
# Ω = PathFinder.NonOscillatoryRegion(RatPhase, ω; Cball = 2π, δball = 1e-3,  Nrays = 16)
# Pexit = PathFinder.exitpoints(RatPhase, Ω)
# CG, CtoG, NodesDict, EdgesList = PathFinder.ContourGraph(RatPhase, a, b, Ω; δODE=0.1, δcoarse=0.01)
valnsd, figs = NumericalSteepestDescent.integrate([-2,-im,2],z->1.0,RatPhase,ω; plot_sd = true)
figs[1]

"""
    Compare Gauss Laguerre convergence
"""

# exit points with ω = 200
ω = 200.
ComplexF64[-1.1147154190894493 + 0.11471541908944906im, -0.8852845809105508 - 0.11471541908944916im, 0.8852845809105508 - 0.11471541908944917im, 1.1147154190894493 + 0.11471541908944916im]

# exit points with ω = 20
ω = 20
exits = ComplexF64[-1.380771508357747 + 0.38077150835774626im, -0.6192284916422535 - 0.38077150835774654im, 0.6192284916422535 - 0.3807715083577466im, 1.3807715083577465 + 0.38077150835774654im]

# exit points with ω = 2
ω = 2
exits = ComplexF64[-1.5639424784471474 + 0.563942478447147im, -0.4360575215528526 - 0.5639424784471472im, 0.4360575215528528 - 0.5639424784471473im, 1.5639424784471472 + 0.5639424784471472im]

g(z) = z + 1/z
dg(z) = 1 - 1/z^2

# Compute Finite contours between exit points and endpoints

# Add SD contour contribution
function h(u,η)
    inv₊(s) = 0.5*(s + sqrt(s^2 - 4))
    inv₋(s) = 0.5*(s - sqrt(s^2 - 4))
    if real(η) > 0 && abs(η) > 1 return inv₊(g(η + im*u)) end
    if real(η) > 0 && abs(η) < 1 return inv₋(g(η + im*u)) end
    if real(η) < 0 && abs(η) > 1 return inv₋(g(η + im*u)) end
    if real(η) < 0 && abs(η) < 1 return inv₊(g(η + im*u)) end
end

nquad = [range(1,20); 100]
sdvals = zeros(ComplexF64, length(exits),length(nquad))
for n in nquad
    x,w = gausslaguerre(n)
    for (i,η) in enumerate(exits)
        # println("η = $η")
        hη = [h(p/ω,η) for p in x]
        dhη = [im / dg(hp) for hp in hη]
        val = sum( w .* dhη ) * cis(ω*g(η))/ω
        # println("val = $val")
        sdvals[i, findfirst(==(n), nquad)] = val
    end
end
refvals = sdvals[:,end]
sderrors = zeros(size(sdvals))
for i in eachindex(exits)
    sderrors[i,:] = abs.(sdvals[i,:] .- refvals[i])
end

fig = Figure()
ax = Axis(fig[1,1], title = "SD contribution convergence", 
        xlabel = "nquad", ylabel = "error", xscale = sqrt, yscale = log10)
for i in eachindex(exits)
    scatterlines!(ax, nquad[1:end-1], sderrors[i,1:end-1], label = "η = $(round(exits[i], digits = 2))")
end
axislegend(ax, position = :lb)
limits!(ax, 1, 20, 1e-14, 1)
fig