"""
    Consider a generic example
    Show quasiSD contour changes with frequency
"""

using PathFinder
using CairoMakie

framerate = 20
Phase = PolynomialPhaseFunction([3,5,6,2,9,5,1,4,1,3])
freqvals = range(0.001,20, step = 1/framerate)
#10 .^ range(-4, 0, step=1/framerate)

frame_iteration(ω) = PathFinder.quasiSDdeformation!(fig, ax, -1,1, Phase, ω;
                                                    umax = 80,
                                                    color_lim = 100,
                                                    resolution = 200,
                                                    set = 1.5)

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -2:1:2, yticks = -2:1:2)
frame_iteration(10000)
fig

record(fig, "generic.gif", freqvals;
        framerate = framerate) do ω
    empty!(ax)
    ax.title = "ω = $(round(ω, digits = 0))"
    frame_iteration(ω)
end

