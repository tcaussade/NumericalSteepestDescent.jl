using PathFinder

ω = 30
f = z -> 1.0
a,b = (1.0, -1/sqrt(2))
a,b = (1.0, -0.8947368421052632)

a,b = 1.0, -0.7894736842105263
a,b = (1.0, -0.9698492462311558)
a,b = 0.001, -1.0

a,b = 1, 0
SqrtPhase = SquareRootPhase(a, b)

val0, figs = nsd(0.0, 1.0,f,SqrtPhase,ω; quadtype = :gaussian, N= 15
,plot_sd=true, plot_graph=true)

figs[1]
figs[2]

ζ = z -> cis(ω * (sqrt(z^2+a^2) + b*z))
abs.(quadgk(ζ, 0.0, 1.0, atol = 1e-14)[1] - val0)


fig = PathFinder.Figure()
ax  = PathFinder.Axis(fig[1, 1], title = "", aspect = PathFinder.DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -5:1:5, yticks = -5:1:5)
PathFinder.quasiSDdeformation!(fig, ax, 0,1, SqrtPhase, 50.;
                                                    umax = 50,
                                                    color_lim = 2)
PathFinder.limits!(-2,2,-2,2)
fig