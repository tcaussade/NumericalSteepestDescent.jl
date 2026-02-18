using PathFinder

"""
    Replicating Fig.3 and Fig.4 in https://p-lpi.github.io
There, the phase is g(z) = 0.5*(x-μ)^2 + α/(1+x^2)

    We are also able to replicate Fig.7 in https://arxiv.org/pdf/1909.04632
"""

ν = 100
α = 1
μ = 0

acoefs = [0.5*μ^2, -μ, 0.5]
ps = [im, -im]
pcoefs = 0.5im*α * [[-1], [1]]

PlasmaLensPhase = RationalPhaseFunction(acoefs, ps, pcoefs)
val, figs = PathFinder.integrate(π,0,z->1.0,PlasmaLensPhase,ν; 
                                infcontour = [true, true],
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]

μvals = range(-0.3,0.3, 50)
intensity = zeros(length(μvals))
for (i,μ) in enumerate(μvals)
    acoefs = [0.5*μ^2, -μ, 0.5]
    PlasmaLensPhase = RationalPhaseFunction(acoefs, ps, pcoefs)
    ψ = integrate(π,0,z->1.0,PlasmaLensPhase,ν; infcontour = [true, true])[1]
    intensity[i] = abs( ψ ) *2 
end

fig = PathFinder.Figure(size = (800,400))
ax = PathFinder.Axis(fig[1, 1], title = "Intensity pattern of a localised lens",  
                xlabel = "μ", ylabel = "Intensity",
                limits = (-0.5,0.5,0,8))
PathFinder.lines!(ax, μvals, intensity)
fig