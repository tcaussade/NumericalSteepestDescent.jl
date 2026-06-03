using NumericalSteepestDescent
using CairoMakie

# Show deformed contour for a given instance of cuspoid integral
framerate = 10
vals = range(-8,8, step = 1/framerate)
# x = 0
y = -1
G(x,y) = PolynomialPhase([0, x, y, 0, 1])

a,b = π/1, 0.0
frame_iteration(x) = NumericalSteepestDescent.quasiSDdeformation!(fig, ax, [a,b], G(x,y), 1;
                                                    umax = 100, color_lim = 100,
                                                    infcontour = [true,true])

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
frame_iteration(first(vals))
limits!(-3,3,-3,3)
fig

record(fig, "animations/cuspoid_x.gif", vals;
        framerate = framerate) do y
    empty!(ax)
    ax.title = "x = $(round(y, digits = 2))"
    frame_iteration(y)
    limits!(-3,3,-3,3)
end

# Show deformed contour for a given instance of catastrophe integral
framerate = 10
xvals = range(-8,8, step = 1/framerate)
y = 1
z = 0
G(x,y,z) = RationalPhase([0,0,(z^2+x),0,2z,0,1],[0.],[[0.0, y^2/12]])

frame_iteration(x) = quasiSDdeformation!(fig, ax, -7π/12,π/12, G(x,y,z), 1;
                                                    umax = 100, color_lim = 10,
                                                    infcontour = [true,true])

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
frame_iteration(first(xvals))
limits!(-2,2,-2,2)
fig

record(fig, "catastrophe.gif", xvals;
        framerate = framerate) do x
    empty!(ax)
    ax.title = "x = $(round(x, digits = 2))"
    frame_iteration(x)
    limits!(-2,2,-2,2)
end