
"""
    Create graph where 
        Vertex are: endpoints, stationary points, entrance/exit points, valleys
        Edges connect two vertices if these are connected by some SD contour outside
        the non-oscillatory region, or Finite contour inside
"""

function ContourGraph(G::AbstractPhaseFunction, a, b, Ω :: Vector{NonOscillatoryBall})
    
    # global rstar = rvalley(G) # compute this only once

    Pendp = [a,b]
    Pexit = exitpoints(G,Ω)
    Pstat = get_Pstat(Ω)
   
    SDpts = ComplexF64.([filter_endpts(a,b,Ω); Pexit])
    # get entrance and valley points
    # Pentr, Dict_entrance, Valleys, Dict_valleys  = tracing_contours(G, SDpts, Ω) 
    γ_to_entrance, γ_to_valley = tracing_contours(G, SDpts, Ω) 

    Valleys = unique([to(γ) for γ in γ_to_valley]) # remove repeated valley if more than one contour goes there
    Pentr = [to(γ) for γ in γ_to_entrance]

    all_nodes = [z for z in [Pexit; [a,b]; Pentr; Valleys; Pstat]] 

    plane_to_graph = Dict{ComplexF64, Int16}() # map complex plane points to graph vertices
    for (i, z) in enumerate(all_nodes)
        plane_to_graph[z] = i
    end

    nvertices = length(plane_to_graph)
    ContourGraph = SimpleGraph(nvertices)
    EdgesList = Dict{Tuple{Int,Int}, ComplexContour}()

    # Part 1: create edges between nodes inside the same ball.
    pts = ComplexF64.([Pexit; Pendp; Pentr; Pstat])
    connect_inside_Ω!(pts, Ω, ContourGraph, plane_to_graph, EdgesList)

    # Part 2: create edges between stationary points with overlapping balls
    connect_stationarypoints!(Ω, ContourGraph, plane_to_graph, EdgesList)

    # Part 3: create edges between exit points and entrance points / valley points
    connect_ball_to_valleyyorentrance!(γ_to_entrance, ContourGraph, plane_to_graph, EdgesList)
    connect_ball_to_valleyyorentrance!(γ_to_valley, ContourGraph, plane_to_graph, EdgesList)
    
    MetaDict = Dict{Symbol, Vector{ComplexF64}}() # useful to plot the graph
    MetaDict[:valleys]   = Valleys
    MetaDict[:entrances] = Pentr
    MetaDict[:exits]     = Pexit
    MetaDict[:statpoint] = Pstat
    MetaDict[:endpoint]  = Pendp

    return ContourGraph, plane_to_graph, MetaDict, EdgesList
end

""" Deal with endpoints 
"""

function filter_endpts(a,b, Ω)
    # a,b are (finite) endpoints
    ComplexF64.([isinΩ(Ω,a) ? [] : a; isinΩ(Ω,b) ? [] : b])
end

""" create connections inside graph """

function connect_inside_Ω!(pts::Vector{ComplexF64}, Ω::Vector{NonOscillatoryBall}, 
                            CG ::SimpleGraph, CtoG::Dict, GtoContour::Dict)
    # create edges between nodes inside the same ball.
    # CtoG : maps points in C to vertex in graph
    # GtoContour : associates graph edge to ComplexContour
    # pts = [Pexit; Pendp; Pentr; Pstat]
    # CG is ContourGraph

    for Ball in Ω
        c, r = centre_and_radius(Ball)
        inside_pts = [z for z in pts if abs(z - c) ≤ r + 1e-14] 
        # added small perturbation here, sometimes the points at the boundary
        # are placed slightly outside the non-oscillatory region.

        for z1 in inside_pts
            for z2 in inside_pts
                if z1 != z2 
                    i1, i2 = CtoG[z1], CtoG[z2]
                    add_edge!(CG, i1,i2)
                    GtoContour[(i1,i2)] = ComplexContour(:finite, z1, z2)
                end
            end
        end
    end

    return
end

function connect_stationarypoints!(Ω::Vector{NonOscillatoryBall}, CG ::SimpleGraph, 
        CtoG::Dict, GtoContour::Dict)
    # create edges between stationary points if their respective balls overlap
    # CG is ContourGraph

    for Ball1 in Ω
        c1,r1 = centre_and_radius(Ball1)
        for Ball2 in Ω
            c2,r2 = centre_and_radius(Ball2)
            if Ball1 != Ball2 && (abs(c1-c2) < r1 +r2)
                i1,i2 = CtoG[c1], CtoG[c2]
                add_edge!(CG, i1,i2)
                GtoContour[(i1,i2)] = ComplexContour(:finite, c1, c2)
            end
        end
    end
    return
end

# function connect_ball_to_valleyyorentrance!(pts::Vector{ComplexF64}, ηdict::Dict, CG ::SimpleGraph, dict::Dict)
#     # create edges between exit points - endpoints and entrance-valleys.
#     # Dict contains η ∈ [Pexit, Pendp] and where it goes
#     for η in pts
#         if haskey(ηdict, η)
#             goesto = ηdict[η]
#             add_edge!(CG, dict[η], dict[goesto])
#         end
#     end
#     return
# end

function connect_ball_to_valleyyorentrance!(γvec, CG ::SimpleGraph, CtoG::Dict, GtoContour::Dict)
    # create edges between exit points - endpoints and entrance-valleys.
    for γ in γvec
        z1, z2 = at(γ), to(γ)
        i1, i2 = CtoG[z1], CtoG[z2]
        add_edge!(CG, i1, i2)
        GtoContour[(i1,i2)] = γ
    end
    return
end

"""
    Plot graph
"""
function plot_ContourGraph(graph::SimpleGraph, Ω::Vector, z_to_G::Dict, metadict :: Dict)
    
    # Place graph nodes in the complex plane
    list = Vector(undef, length(z_to_G))

    # Add colors to nodes
    # Create a "Vector of colors" to pass an input to graph
    colors = Vector{Any}(undef, length(z_to_G))
    for z in metadict[:statpoint]
        i = z_to_G[z]
        colors[i] = colorant"red"
        list[i]   = reim(z)
    end
    for z in metadict[:endpoint]
        i = z_to_G[z]
        colors[i] = colorant"orange"
        list[i]   = reim(z)
    end
    for z in metadict[:exits]
        i = z_to_G[z]
        colors[i] = colorant"purple"
        list[i]   = reim(z)
    end
    for z in metadict[:valleys]
        i = z_to_G[z]
        colors[i] = colorant"blue"
        list[i]   = reim(2 * z)
    end
    for z in metadict[:entrances]
        i = z_to_G[z]
        colors[i] = colorant"green"
        list[i]   = reim(z)
    end


    fig,ax,p = graphplot(graph,
        node_color = colors, node_size = 20)

    mylayout(_) = list
    p.layout = mylayout

    # add non-oscillatory balls
    for Ball in Ω
        c,r = centre_and_radius(Ball)
        arc!(ax, Point2f(reim(c)), r, 0, 2π, 
            color = :gray)
    end

    hidespines!(ax)
    hidedecorations!(ax)
    ax.aspect = DataAspect()
    ax.backgroundcolor

    return fig

end