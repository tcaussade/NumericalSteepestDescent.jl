
# Plot SD lansdscape and lines of constant real part equal to stationary points.
function plot_landscape(G :: AbstractPhase; ran = 10, color_lim = 100)
    nx = 200
    g(z) = evalphase(z, G)
    ξ = stationary_points(G)
    x = range(-ran,ran, length = nx)
    y = range(-ran,ran, length = nx)
    X = [x for x in x for _ in y]
    Y = [y for _ in x for y in y]
    Z = [g(x+im*y) for x in x for y in y]

    rlevels = unique(real.(g.(ξ)))

    # color_lim = 300 # maximum(imag(Z)) / 2
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Re(z)", ylabel = "Im(z)", title = "SD landscape", aspect = 1)
    levelset = contourf!(ax,X,Y,-imag.(Z); levels = range(-color_lim, color_lim, 20), 
                         colormap = :jet, extendlow = :auto, extendhigh = :auto)
    contour!(ax,X,Y,real.(Z); levels = rlevels, color = :black, linewidth = 1, linestyle = :dash)
    scatter!(ax, reim.(ξ), color = :red, markersize = 12, marker = :star5) 
    # scatter!(ax, reim.([η]), color = :blue, markersize = 12, marker = :circle)
    Colorbar(fig[1,2], levelset)
    return fig
end