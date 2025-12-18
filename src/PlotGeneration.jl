"""
    Functionalities to plot the quasi-SD contour deformation
"""

function plot_SDcontours(G::AbstractPhaseFunction, γ::Vector{ComplexContour}, Ω)

    @show γ
    resolution = 200

    xmin = -2
    xmax = +2
    ymin = -2
    ymax = +2
    x = range(xmin,xmax, resolution)
    y = range(xmin,xmax, resolution)
    θ = range(0, 2π, resolution)

    u = collect(range(0,20,resolution)) # used for SD contours  
    # t = collect(range(-1,1,2))          # used for finite contours


    fig = Figure()
    ax = Axis(fig[1, 1], title = "Quasi-SD deformation", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im")
              #xticks = -xmin:1:xmax, yticks = -xmin:1:)
    
    # plot levelset of phase function
    X = [x for x in x for _ in y]
    Y = [y for _ in x for y in y]
    Z = [evalphase(G, x+im*y) for x in x for y in y]
    levelset = contourf!(ax,X,Y,-imag.(Z); levels = 20, colormap = :balance)
    # contour!(ax,X,Y,real.(Z); levels = 11, color = :black, linewidth = 1, linestyle = :dash)
        
    # plot non-oscillatory region
    for ball in Ω # Display Non-oscillatory region(s)
        c,r = centre_and_radius(ball)
        zb = c .+ r*cispi.(θ)
        lines!(ax,reim.(zb), color = :gray)
    end

    # add contours
    for c in γ
        if contour_type(c) == :infiniteSD
            hη = points_on_SDcontour(at(c), G.p, G.dp, u; δfine = 1e-2)
            lines!(ax, reim.(hη); color = :blue, linewidth = 2)
        elseif contour_type(c) == :finite
            lines!(reim.([at(c), to(c)]); color = :red, linewidth = 2)
        elseif contour_type(c) == :finiteSD
            U = im*(G.p(at(c)) - G.p(to(c)))
            u_tmp = u * U/20
            hη = points_on_SDcontour(at(c), G.p, G.dp, u_tmp; δfine = 1e-2)
            lines!(ax, reim.(hη); color = :green, linewidth = 2)
        end
    end

    # add stationary points and endpoints
    scatter!(ax, reim.(G.ξ), color = :red)
    scatter!(ax, reim.([at(γ[1]), at(γ[end])]), color = "black")
    limits!(xmin,xmax,ymin,ymax)
    Colorbar(fig[1,2], levelset)
    return fig
end
