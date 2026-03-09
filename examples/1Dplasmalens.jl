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

"""
    Replicating Fig.2 in "Oscillatory path integrals for radio astronomy"
    (Feldbrugge, Pen, Turol, 2023)

    The phase is g(z) = 0.5*(x-μ)^2 + 0.5*α/(1+x^2) where (see eq.(1))
        - μ is position
        - α is lens strength
"""

μvals = range(-0.6,0.6, 31)
αvals = range(1e-5,4,30)
ν = 10

I = zeros(length(μvals),length(αvals))
ps = [im, -im]
for (j,α) in enumerate(αvals)
    pcoefs = 0.25im*α * [[-1], [1]]
    for (m,μ) in enumerate(μvals)
        println("Evaluating at (μ,α) = ($μ,$α)")
        acoefs = [0.5*μ^2, -μ, 0.5]
        G = RationalPhaseFunction(acoefs, ps, pcoefs)
        ψ = integrate(π,0,z->1.0,G,ν; infcontour = [true, true])
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
fig