"""
    High-order coalescent phase quasi SD contour deformation
    g(z) = 1/7 *z^7 - r^6 * z
    with r>0 
"""

using PathFinder
using CairoMakie

framerate = 100

P = 6
G(r) = PolynomialPhaseFunction([0.0; -r^P; zeros(P-1); 1.0/(P+1)])
rvals = range(1e-4, 1e0, step = 1/framerate)
#10 .^ range(-4, 0, step=1/framerate)

frame_iteration(r) = PathFinder.quasiSDdeformation!(fig, ax, -1,1, G(r), 1000.0;
                                                    umax = 10,
                                                    color_lim = 8,
                                                    resolution = 500,
                                                    set = 1.5)

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -2:0.5:2, yticks = -2:0.5:2)
frame_iteration(first(rvals))
fig

record(fig, "coalescence.gif", rvals;
        framerate = framerate) do r
    empty!(ax)
    ax.title = "r = $(round(r, digits = 1))"
    frame_iteration(r)
end

