using PathFinder

ω = 30
f = z -> 1.0
a,b = (1.0, -1/sqrt(2))
a,b = (1.0, -0.8947368421052632)

a,b = 1.0, -0.7894736842105263
a,b = (1.0, -0.9698492462311558)
a,b = 1e-6, -1.0

SqrtPhase = SquareRootPhaseFunction(a, b)

val0, figs = integrate(0.0, 1.0,f,SqrtPhase,ω; quadtype = :gaussian, N=15,
                plot_sd=true, plot_graph=true)

figs[1]
figs[2]

ζ = z -> cis(ω * (sqrt(z^2+a^2) + b*z))
quadgk(ζ, 0.0, 1.0, atol = 1e-14)[1]