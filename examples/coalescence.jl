"""
    High-order coalescent phase quasi SD contour deformation
    g(z) = 1/7 *z^7 - r^6 * z
    with r>0 
"""

using PathFinder
using CairoMakie

framerate = 20

P = 6
G(r) = PolynomialPhaseFunction([0.0; -r^P; zeros(P-1); 1.0/(P+1)])
rvals = 10 .^ range(-4, 0, length=framerate)

frame_iteration(r) = PathFinder.quasiSDdeformation!(fig, ax, -1,1, G(r), 1000.0;
                                                    umax = 150,
                                                    color_lim = 8)

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -2:0.5:2, yticks = -2:0.5:2)
frame_iteration(first(rvals))
limits!(-1.5,1.5,-1.5,1.5)
fig

record(fig, "coalescence.gif", rvals;
        framerate = framerate) do r
    empty!(ax)
    ax.title = "r = $(round(r, digits = 1))"
    frame_iteration(r)
    limits!(-1.5,1.5,-1.5,1.5)
end
