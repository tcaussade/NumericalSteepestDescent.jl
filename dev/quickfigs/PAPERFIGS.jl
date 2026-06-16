using NumericalSteepestDescent
using CairoMakie

ω = 5

Mon = PolynomialPhase([0,0,1])
i1,fig = nsd([cis(-π/4),1],x -> 1,Mon,ω; 
                N = 50, infcontour = [false, false], 
                plot_sd = true)
i1
limits!(-2,2,-2,2)
fig[1]

a=-0.0
Airy = PolynomialPhase([0,-3a,0,1])
ω = 5
γ0, infcontour = [5π/6, π/6], [true,true]
fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
quasiSDdeformation!(fig, ax, γ0, Airy ,ω; infcontour, umax = 150,
                                                     color_lim = 100)
limits!(-4,4,-4,4)
fig