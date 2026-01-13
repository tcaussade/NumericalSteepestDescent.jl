"""
    Functionalities to plot the quasi-SD contour deformation
"""

function plot_SDcontours(G::AbstractPhaseFunction, γ::Vector{ComplexContour}, Ω, γall::Vector{ComplexContour};
        infcontour, inftol,
        color_lim = 100)

        fig = Figure()
        ax = Axis(fig[1, 1], title = "Quasi-SD deformation", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im")
              #xticks = -xmin:1:xmax, yticks = -xmin:1:)
        
        plot_SDcontours!(fig, ax, G,γ,Ω, γall; infcontour, inftol, color_lim)
end

function plot_SDcontours!(fig, ax, G::AbstractPhaseFunction, γ::Vector{ComplexContour}, Ω, γall::Vector{ComplexContour};
        infcontour, inftol,
        umax = 50, # control how far tracing a contour for plots
        color_lim  = 300, # control color limits of the colorbar
        resolution = 200, # increasing this paramter improves image quality
        set        = 10 # plotsize
        )

    xmin = -set
    xmax = +set
    ymin = -set
    ymax = +set
    x = range(xmin,xmax, resolution)
    y = range(xmin,xmax, resolution)
    θ = range(0, 2, resolution)

    # umax = 50
    u = collect(range(0,umax,10*resolution)) # used for SD contours  

    # plot levelset of phase function
    X = [x for x in x for _ in y]
    Y = [y for _ in x for y in y]
    Z = [evalphase(G, x+im*y) for x in x for y in y]
    # color_lim = 300 # maximum(imag(Z)) / 2
    levelset = contourf!(ax,X,Y,-imag.(Z); levels = range(-color_lim, color_lim, 20), 
                         colormap = :balance, extendlow = :auto, extendhigh = :auto)
    # contour!(ax,X,Y,real.(Z); levels = 11, color = :black, linewidth = 1, linestyle = :dash)
        
    # plot non-oscillatory region
    for ball in Ω # Display Non-oscillatory region(s)
        c,r = centre_and_radius(ball)
        zb = c .+ r*cispi.(θ)
        lines!(ax,reim.(zb), color = :gray)
    end

    # add contours of the quasi-SD deformation
    g(z)  = evalphase(G,z)
    dg(z) = evalphase_derivative(G,z)
    for c in γall
        if abs(at(c)) > inftol 
            @show abs(at(c))
            continue 
        end
        lw =  c in γ ? 3 : 1 # use wider line for SD contours on shortest path
        if contour_type(c) == :infiniteSD
            hη = points_on_SDcontour(at(c), G, u; δfine = 1e-12)
            lines!(ax, reim.(hη); color = :blue, linewidth = lw)
        elseif contour_type(c) == :finiteSD
            U = im*(G.p(at(c)) - G.p(to(c)))
            u_tmp = u * U/umax
            hη = points_on_SDcontour(at(c), G, u_tmp; δfine = 1e-6)
            lines!(ax, reim.(hη); color = :green, linewidth = lw)
        end
    end
    for c in γ
        if contour_type(c) == :finite
            lines!(reim.([at(c), to(c)]); color = :red, linewidth = 3)
        end
    end

    # add stationary points and endpoints
    scatter!(ax, reim.(G.ξ), color = :red)
    !infcontour[1] ? scatter!(ax, reim.([at(γ[1])]), color = "black") : nothing
    display = contour_type(γ[end]) == :finite ? to(γ[end]) : at(γ[end])
    !infcontour[2] ? scatter!(ax, reim.([display]), color = "black") : nothing
    # scatter!(ax, reim.([at(γ[1]), at(γ[end])]), color = "black")
    limits!(xmin,xmax,ymin,ymax)
    Colorbar(fig[1,2], levelset)
    return fig
end
