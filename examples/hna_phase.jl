
"""
    Phase encountered in a BEM matrix for HNA methods
"""

using NumericalSteepestDescent
using CairoMakie

a = 0.1

framerate = 60

bvals = range(-0.99,0.99, step = 1/framerate)
G(a,b) = SquareRootPhaseFunction(a,b)

frame_iteration(b) = PathFinder.quasiSDdeformation!(fig, ax, [0,1], G(a,b), 100;
                                                    umax = 20,
                                                    color_lim = 2)

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
frame_iteration(first(bvals))
limits!(-5,5,-5,5)

record(fig, "hnaphase_"*"$a"*".gif", bvals;
        framerate = framerate) do b
    empty!(ax)
    ax.title = "b = $(round(b, digits = 2))"
    frame_iteration(b)
    limits!(-5,5,-5,5)
end
