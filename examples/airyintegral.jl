
"""
    Airy phase quasi SD contour deformation
"""

using NumericalSteepestDescent
using CairoMakie

framerate = 20

xvals = range(-5,5, step = 1/framerate)
G(x) = PolynomialPhase(-im*[0,-x,0,1/3])

frame_iteration(x) = PathFinder.quasiSDdeformation!(fig, ax, [-π/3,π/3], G(x), 1.0; infcontour = [true, true],
                                                    umax = 150,
                                                    color_lim = 100)

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
frame_iteration(first(xvals))
limits!(-6,6,-6,6)

record(fig, "airy.gif", xvals;
        framerate = framerate) do x
    empty!(ax)
    ax.title = "x = $(round(x, digits = 1))"
    frame_iteration(x)
    limits!(-6,6,-6,6)
end
