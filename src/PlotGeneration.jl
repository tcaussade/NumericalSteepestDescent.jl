"""
    Functionalities to plot the quasi-SD contour deformation
"""

function plot_quasiSDdeformation(G::AbstractPhaseFunction, γ::Vector{ComplexContour}, Ω)

    fig = Figure()
    ax = Axis(fig[1, 1], title = "Quasi-SD deformation", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -3:1:3, yticks = -3:1:3)

    plot_quasiSDdeformation(fig, G, γ, Ω)
end


function plot_quasiSDdeformation(fig, G::AbstractPhaseFunction, γ::Vector{ComplexContour}, Ω)

    u = range(0,15,400)
    t = range(-1,1,50)
    
    for ball in Ω # Display Non-oscillatory region(s)
        zb = ball[1] .+ ball[2]*cispi.(t)
        lines!(reim.(zb), color = :gray)
    end

    for c in γ
        if contour_type(c) == :infiniteSD
            lines!(reim.(c.parametrisation.(u)); color = :blue, linewidth = 2)
        elseif contour_type(c) == :finite
            lines!(reim.(c.parametrisation.(t)); color = :red, linewidth = 2)
        elseif contour_type(c) == :finiteSD
            lines!(reim.(c.parametrisation.(t)); color = :green, linewidth = 2)
        end
    end
    scatter!(reim.(G.ξ), color = :red)
    scatter!(reim.([0.0im,1.0+0im]), color = "black")
    limits!(-2,2,-2,2)
    return fig
end