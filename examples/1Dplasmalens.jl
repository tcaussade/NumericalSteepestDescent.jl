using NumericalSteepestDescent
using CairoMakie

"""
    Replicating Fig.3 and Fig.4 in https://p-lpi.github.io
There, the phase is g(z) = 0.5*(x-μ)^2 + α/(1+x^2)

    We are also able to replicate Fig.7 in https://arxiv.org/pdf/1909.04632
"""


ν = 40
α = 2
μ = 0

acoefs = [0.5*μ^2, -μ, 0.5]
ps = [im, -im]
pcoefs = 0.25im*α * [[-1], [1]]

PlasmaLensPhase = RationalPhaseFunction(acoefs, ps, pcoefs)
val, figs = PathFinder.integrate([π,0],z->1.0,PlasmaLensPhase,ν; 
                                infcontour = [true, true],
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]

""" Animation for a range of parameters """
framerate = 80
μvals = range(-0.6,0.6, step = 1/framerate)
α = 5

PlasmaPhase(μ) = RationalPhaseFunction([0.5*μ^2, -μ, 0.5], [im, -im], 0.25im*α * [[-1], [1]])
frame_iteration(μ) = PathFinder.quasiSDdeformation!(fig, ax, [π,0.0], PlasmaPhase(μ), ν; 
                                                    infcontour = [true, true],
                                                    umax = 50,
                                                    color_lim = 10)
fig = PathFinder.Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
frame_iteration(first(μvals))
limits!(-4,4,-4,4)
fig

record(fig, "animations/plasma.gif", μvals;
        framerate = framerate) do μ
    empty!(ax)
    ax.title = "α = $α, μ = $(round(μ, digits = 1))"
    frame_iteration(μ)
    limits!(-4,4,-4,4)
end


"""
    Replicating Fig.2 in "Oscillatory path integrals for radio astronomy"
    (Feldbrugge, Pen, Turol, 2023)

    The phase is g(z) = 0.5*(x-μ)^2 + 0.5*α/(1+x^2) where (see eq.(1))
        - μ is position
        - α is lens strength
"""

μvals = range(-0.6,0.6, 51)
αvals = range(1e-5,6,51)
ν = 20

I = zeros(length(μvals),length(αvals))
ps = [im, -im]
for (j,α) in enumerate(αvals)
    pcoefs = 0.25im*α * [[-1], [1]]
    for (m,μ) in enumerate(μvals)
        # println("Evaluating at (μ,α) = ($μ,$α)")
        acoefs = [0.5*μ^2, -μ, 0.5]
        G = RationalPhaseFunction(acoefs, ps, pcoefs)
        ψ = integrate([π,0],z->1.0,G,ν; infcontour = [true, true])
        I[m,j] = abs(ψ * sqrt(-im*ν/(2π)))^2
    end
end

fig = PathFinder.Figure()
ax = PathFinder.Axis(fig[1, 1], title = "Intensity pattern of a localised lens (ν=$ν)",  
                xlabel = "μ", ylabel = "α")
levelset = PathFinder.contourf!(ax,μvals,αvals,I; levels = range(0, 5, 100),
                    colormap = :jet, # colormap = :hot
                    extendlow = :auto, extendhigh = :auto)
PathFinder.Colorbar(fig[1,2], levelset)
limits!(first(μvals), last(μvals), first(αvals), last(αvals))
fig