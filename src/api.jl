"""
    Generate QuasiSDcontour given two finite endpoints
"""


function QuasiSDcontour(G::AbstractPhaseFunction, a, b, Ω :: Vector{NonOscillatoryBall})
    # a,b are (finite) endpoints
    
    CG, dict, nodes = ContourGraph(G, a, b, Ω)
    
    sd = dijkstra_shortest_paths(graph, dict[a])

    sd.dists

    
    return 
end


"""
    Evaluate over a given quasi-SD contour
"""

function integrate_nsp(f, G::AbstractPhaseFunction, γ::Vector{ComplexContour}, ω, N)
    # evaluate integral along quasi-SD contour deformation
    @info "Using N=$N quadrature points for all contour segments. \nThese are being recomputed at each call!"
    x1,w1 = gausslegendre(N)
    x2,w2 = gausslaguerre(N)
    x3,w3 = x1,w1 # gausslegendre(N);

    I = 0.0 + 0im
    for c in γ
        ± = i -> contour_orientation(c) == :positive ? +(i) : -(i)
        if contour_type(c) == :finite
            I += ± eval_finite(f, G, c, ω, x1, w1)
        elseif contour_type(c) == :infiniteSD
            I += ± eval_infiniteSDpath(f, G, c, ω, x2, w2)
        elseif contour_type(c) == :finiteSD
            I += ± eval_finiteSD(f, G, c, ω, x3, w3)
        end
    end
    return I
end
