"""
    Compute the oscillatory integral with
        - finite endpoints (a,b), or infinite endpoints with angle (a,b) if specified
        - amplitude function f
        - phase function G
        - frequency ω
    Returns the value and a vector of figures
"""

function integrate(a, b, f::Function, G::AbstractPhaseFunction, ω; 
        # default parameters
        Cball = 2π,    # control maximum number of oscillations on non-oscillatory bals
        Nrays = 16,    # number of rays used to determine ball radius
        δball = 1e-3,  # determine when overlapping balls should be amalgamated
        δODE  = 0.1,   # local step size in ODE solver for SD path tracking
        δcoarse= 0.01, # corrector tolerance in SD tracking
        δfine = 1e-13, # tolerance to compute weights and nodes along SD contours
        δquad = 1e-16, # used for truncation and to determine when a contour should be dropped

        infcontour = [false, false], # specify if endpoints are at infinity

        # quadrature
        quadtype = :gaussian, # specify quadrature type [:gaussian (default), :adaptive]   
        N = 25,        # number of quadrature points (ignored if quadtype == :adaptive)
        atol = 1e-10, # tolerance for absolute error in adaptive quadrature. 
        # WARNING: setting rtol is unreliable for integrals with small values 

        # produce plots 
        plot_graph = false, # if true, returns the graph plot
        plot_sd    = false, # if true, plots the chosen quasi-SD contour for evaluation
        inftol = 1e7,  # used to discard tracing from points that are too large in plots
        )
  
    # place endpoints at infinity if specified
    a = infcontour[1] ? endpoint_at_valley(a, G) : a
    b = infcontour[2] ? endpoint_at_valley(b, G) : b
    
    Ω = NonOscillatoryRegion(G, ω; Cball, δball,  Nrays)

    CG, CtoG, NodesDict, EdgesList = ContourGraph(G, a, b, Ω; δODE, δcoarse)
    # CG is ContourGraph, CtoG maps complex plane points to graph vertices
    # NodesDict contains the different types of nodes in the graph
    # EdgesList maps graph edges to ComplexContours
    a,b = NodesDict[:endpoint]
    sd_edges = a_star(CG, CtoG[a], CtoG[b]) # find shortest path

    # choose quadrature and precompute weights and nodes
    if quadtype == :gaussian
        qleg = gausslegendre(N)
        qlag = gausslaguerre(N)
        if G isa SquareRootPhaseFunction 
            qsin = gausslegendre(Int(floor(N/layersnumber(G.a))) + 1) 
        else qsin = nothing end
        quad = (qleg = qleg, qlag = qlag, qsingular = qsin)
    elseif quadtype == :adaptive
    else
        @error "quadtype is not in (:gaussian, :adaptive)"
    end

    # Evaluate each contour on the shortest path
    S = zero(ComplexF64)
    γtot = Vector{ComplexContour}()
    for e in sd_edges 
        i1,i2 = e.src, e.dst
        if haskey(EdgesList, (i1, i2))
            γ = EdgesList[(i1,i2)]; push!(γtot, γ)
            if quadtype == :gaussian
                S += integrate(γ, f, G, ω, quad; δfine, δquad)
            elseif quadtype == :adaptive
                S += integrate(γ, f, G, ω; δfine, δquad, atol)
            end
        else # the contour is traversed in the opposite direction
            γ = EdgesList[(i2,i1)]; push!(γtot, γ)
            if quadtype == :gaussian
                S -= integrate(γ, f, G, ω, quad; δfine, δquad)
            elseif quadtype == :adaptive
                S -= integrate(γ, f, G, ω; δfine, δquad, atol)
            end
        end
    end

    # Add contribution of residues if necessary
    # S += add_residues(G,γtot)

    if abs(S) > 1e12
        @warn "Value of the integral is large"
    end

    if isempty(γtot) 
        @warn "The graph is not connected between endpoints!"
        return nothing, [plot_ContourGraph(CG, Ω, CtoG, NodesDict)]
    end
    
    # for γ in γtot
    #     @show γ
    # end

    if plot_graph || plot_sd
        figs = []
        if plot_graph
            fig1 = plot_ContourGraph(CG, Ω, CtoG, NodesDict)
            push!(figs, fig1)
        end

        if plot_sd 
            γall = Vector{ComplexContour}() # contains all traced contours
            for ηi in [NodesDict[:exits]; NodesDict[:endpoint]]
                for ηj in [NodesDict[:valleys]; NodesDict[:entrances]; NodesDict[:poles]]
                    i = CtoG[ηi]
                    j = CtoG[ηj]
                    if haskey(EdgesList, (i,j)) push!(γall, EdgesList[(i,j)]) end
                end
            end
            fig2 =  plot_SDcontours(G,γtot, Ω, γall; infcontour, inftol)
            push!(figs, fig2)
        end

        return S, figs
    end
    #else, if no plot is needed...
    return S
end


function endpoint_at_valley(θ, G::AbstractPhaseFunction)
    # place endpoint at valley if specified as endpoint at infinity
    v = valleyangle(θ, G)
    if v isa Nothing @warn "endpoint with θ=$(θ/π)π  is not in valley region" end
    # THIS IS A PATCH FIX TO PUT THE POINT OUTSIDE NonOscillatoryRegion!!
    return (rstar_valley(G) + 5) * cis(v) # place far away in valley direction
end


""" 
    Show QuasiSD contour deformation 
    This method is useful for creating gif animations.
"""
function quasiSDdeformation!(fig::Figure,ax::Axis, a,b, G::AbstractPhaseFunction, ω; 
                         infcontour = [false, false], 
                         Cball = 2π,
                         inftol = 1e6,
                         umax = 10,
                         color_lim = 200,
                         resolution = 200,
                         set = 10)

    a = infcontour[1] ? endpoint_at_valley(a,G) : a
    b = infcontour[2] ? endpoint_at_valley(b,G) : b

    Ω = NonOscillatoryRegion(G, ω; Cball, δball=1e-3,  Nrays=16)
    CG, CtoG, NodesDict, EdgesList = ContourGraph(G, a, b, Ω; δODE=0.1, δcoarse=0.01)
    a,b = NodesDict[:endpoint]
    sd_edges = a_star(CG, CtoG[a], CtoG[b]) # find shortest path

    γtot = Vector{ComplexContour}()
    for e in sd_edges i1,i2 = e.src, e.dst
        if haskey(EdgesList, (i1, i2)) γ = EdgesList[(i1,i2)]; push!(γtot, γ)
        else γ = EdgesList[(i2,i1)]; push!(γtot, γ)
        end
    end

    γall = Vector{ComplexContour}() 
    for ηi in [NodesDict[:exits]; NodesDict[:endpoint]]
        for ηj in [NodesDict[:valleys]; NodesDict[:entrances]]
            i,j = CtoG[ηi], CtoG[ηj]
            if haskey(EdgesList, (i,j)) push!(γall, EdgesList[(i,j)]) end
        end
    end

    return plot_SDcontours!(fig, ax, G,γtot, Ω, γall; infcontour, inftol, umax, color_lim, resolution, set)
end

function showContourGraph!(fig::Figure,ax::Axis, a,b, G::AbstractPhaseFunction, ω; 
                         infcontour = [false, false], 
                         Cball = 2π,
                         inftol = 1e6,
                         umax = 10,
                         color_lim = 200,
                         resolution = 200,
                         set = 10)

    a = infcontour[1] ? endpoint_at_valley(a,G) : a
    b = infcontour[2] ? endpoint_at_valley(b,G) : b

    Ω = NonOscillatoryRegion(G, ω; Cball, δball=1e-3,  Nrays=16)
    CG, CtoG, NodesDict, _ = ContourGraph(G, a, b, Ω; δODE=0.1, δcoarse=0.01)
    plot_ContourGraph(CG, Ω, CtoG, NodesDict)
end