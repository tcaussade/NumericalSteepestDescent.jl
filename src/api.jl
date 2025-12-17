"""
    Generate QuasiSDcontour given two finite endpoints
"""

function integrate(a, b, f::Function, G::AbstractPhaseFunction, ω; 
        N = 25, Cball = 2π)
    # a,b are (finite) endpoints
    
    Ω = NonOscillatoryRegion(G, Cball, ω)

    CG, CtoG, MetaDict, EdgesList = ContourGraph(G, a, b, Ω)
    
    sd_edges = a_star(CG, CtoG[a], CtoG[b]) # find shortest path

    # Evaluate each contour on the shortest path
    global xleg, wleg = gausslegendre(N)
    global xlag, wlag = gausslaguerre(N)

    S = zero(ComplexF64)
    for e in sd_edges
        i1,i2 = e.src, e.dst

        if haskey(EdgesList, (i1, i2))
            γ = EdgesList[(i1,i2)]
            x,w = choose_quadrature(γ)
            S += integrate(γ, f, G, ω, x, w)
        else
            γ = EdgesList[(i2,i1)] # the contour is traversed in the opposite direction
            x,w = choose_quadrature(γ)    
            S -= integrate(γ, f, G, ω, x, w)
        end
    end
    return S
end

function choose_quadrature(γ)
    if contour_type(γ) == :infiniteSD # Choose quadrature nodes
        return x,w = xlag,wlag
    else 
        return x,w = xleg,wleg 
    end
end