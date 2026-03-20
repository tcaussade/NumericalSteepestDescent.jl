using NumericalSteepestDescent
using CairoMakie

# Show deformed contour for a given instance of cuspoid integral
framerate = 10
yvals = range(-8,8, step = 1/framerate)
x = 0
G(x,y) = PolynomialPhaseFunction([0, x, y, 0, 1])

frame_iteration(y) = PathFinder.quasiSDdeformation!(fig, ax, [a,b], G(x,y), 1;
                                                    umax = 100, color_lim = 100,
                                                    infcontour = [true,true])

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
frame_iteration(first(xvals))
limits!(-3,3,-3,3)
fig

record(fig, "cuspoid.gif", yvals;
        framerate = framerate) do y
    empty!(ax)
    ax.title = "y = $(round(y, digits = 2))"
    frame_iteration(y)
    limits!(-3,3,-3,3)
end

# Show deformed contour for a given instance of catastrophe integral
framerate = 10
xvals = range(-8,8, step = 1/framerate)
y = 1
z = 0
G(x,y,z) = RationalPhaseFunction([0,0,(z^2+x),0,2z,0,1],[0.],[[0.0, y^2/12]])

frame_iteration(x) = PathFinder.quasiSDdeformation!(fig, ax, -7π/12,π/12, G(x,y,z), 1;
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