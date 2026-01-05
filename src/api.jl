"""
    Generate QuasiSDcontour given two finite endpoints
"""

function integrate(a, b, f::Function, G::AbstractPhaseFunction, ω; 
        # default parameters
        N = 25,        # number of quadrature points
        Cball = 2π,    # control maximum number of oscillations on non-oscillatory bals
        Nrays = 16,    # number of rays used to determine ball radius
        δball = 1e-3,  # determine when overlapping balls should be amalgamated
        δODE  = 0.1,   # local step size in ODE solver for SD path tracking
        δcoarse= 0.01, # corrector tolerance in SD tracking
        δfine = 1e-13, # tolerance to compute weights and nodes along SD contours
        δquad = 1e-16, # used for truncation and to determine when a contour should be dropped

        # produce plots 
        plot_graph = false, # if true, returns the graph plot
        plot_sd    = false, # if true, plots the chosen quasi-SD contour for evaluation
        )
    # a,b are (finite) endpoints
    
    Ω = NonOscillatoryRegion(G, ω; Cball, δball,  Nrays)

    CG, CtoG, NodesDict, EdgesList = ContourGraph(G, a, b, Ω; δODE, δcoarse)
    a,b = NodesDict[:endpoint]
    sd_edges = a_star(CG, CtoG[a], CtoG[b]) # find shortest path

    # Evaluate each contour on the shortest path
    global xleg, wleg = gausslegendre(N)
    global xlag, wlag = gausslaguerre(N)

    S = zero(ComplexF64)
    γtot = Vector{ComplexContour}()
    for e in sd_edges
        i1,i2 = e.src, e.dst
        if haskey(EdgesList, (i1, i2))
            γ = EdgesList[(i1,i2)]
            push!(γtot, γ)
            x,w = choose_quadrature(γ)
            S += integrate(γ, f, G, ω, x, w; δfine, δquad)
        else
            γ = EdgesList[(i2,i1)] # the contour is traversed in the opposite direction
            push!(γtot, γ)
            x,w = choose_quadrature(γ)    
            S -= integrate(γ, f, G, ω, x, w; δfine, δquad)
        end
    end

    fig1 = plot_graph ? plot_ContourGraph(CG, Ω, CtoG, NodesDict) : nothing
    fig2 = plot_sd ? plot_SDcontours(G,γtot, Ω) : nothing
    figs = [fig1, fig2]

    return S, figs
end

function choose_quadrature(γ)
    if contour_type(γ) == :infiniteSD # Choose quadrature nodes
        return x,w = xlag,wlag
    else 
        return x,w = xleg,wleg 
    end
end