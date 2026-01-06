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

        infcontour = [false, false], # specify if endpoints are at infinity
        # quad = :gaussian, # specify quadrature type [:gaussian (default), :adaptive]   

        # produce plots 
        plot_graph = false, # if true, returns the graph plot
        plot_sd    = false, # if true, plots the chosen quasi-SD contour for evaluation
        )
  
    # place endpoints at infinity if specified
    a = infcontour[1] ? endpoint_at_valley!(G, a) : a
    b = infcontour[2] ? endpoint_at_valley!(G, b) : b
    
    Ω = NonOscillatoryRegion(G, ω; Cball, δball,  Nrays)

    CG, CtoG, NodesDict, EdgesList = ContourGraph(G, a, b, Ω; δODE, δcoarse)
    # CG is ContourGraph, CtoG maps complex plane points to graph vertices
    # NodesDict contains the different types of nodes in the graph
    # EdgesList maps graph edges to ComplexContours
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

    # @assert !isempty(γtot) "The graph is not connected between endpoints!"

    γall = Vector{ComplexContour}() # contains all traced contours
    for ηi in NodesDict[:exits] 
        for ηj in [NodesDict[:valleys]; NodesDict[:entrances]]
            i = CtoG[ηi]
            j = CtoG[ηj]
            if haskey(EdgesList, (i,j)) push!(γall, EdgesList[(i,j)]) end
        end
    end
    
    fig1 = plot_graph ? plot_ContourGraph(CG, Ω, CtoG, NodesDict) : nothing
    fig2 = plot_sd ? plot_SDcontours(G,γtot, Ω, γall; infcontour) : nothing
    figs = [fig1, fig2]

    return S, figs
end

function endpoint_at_valley!(G::AbstractPhaseFunction, θ)
    # place endpoint at valley if specified as endpoint at infinity
    v = goes_to_valley(G, θ)
    if v isa Nothing @warn "endpoint with θ=$(θ/π)π  is not in valley region" end
    # THIS IS A PATCH FIX TO PUT THE VALLEY OUTSIDE NonOscillatoryRegion!!
    return (G.rStar + 5) * cis(v) # place far away in valley direction
end

function choose_quadrature(γ)
    if contour_type(γ) == :infiniteSD # Choose quadrature nodes
        return x,w = xlag,wlag
    else 
        return x,w = xleg,wleg 
    end
end