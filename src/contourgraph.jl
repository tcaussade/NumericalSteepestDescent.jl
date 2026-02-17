
"""
    Create graph where 
        Vertex are: endpoints, stationary points, entrance/exit points, valleys
        Edges connect two vertices if these are connected by some SD contour outside
        the non-oscillatory region, or Finite contour inside
"""

function ContourGraph(G::AbstractPhaseFunction, a, b, Ω :: Vector{NonOscillatoryBall};
                        δODE, δcoarse)
    
    # Flag all nodes
    NodesDict = Dict{Symbol, Vector{ComplexF64}}()
    
    # add small perturbation to distinguish between endpoints and stationary points when they coincide
    NodesDict[:endpoint]  = [a,b]        .+ im* eps() * 10
    NodesDict[:statpoint] = get_Pstat(Ω)
    NodesDict[:exits] = exitpoints(G,Ω) #.+ rand()*eps()
    exits_outside_Ω = _filter_exits(G, Ω, NodesDict[:exits]) # used to filter manually placed exit points
    # trace all possible SD contours
    aδ, bδ = NodesDict[:endpoint]
    endpoints_outside_Ω = ComplexF64.([isinΩ(Ω,a) ? [] : aδ; isinΩ(Ω,b) ? [] : bδ])

    endpoints_outside_Ω, exits_outside_Ω
    trace_from = [endpoints_outside_Ω; exits_outside_Ω] 
    γ_to_entrance, γ_to_valley, γ_to_pole = tracing_contours(G, trace_from, Ω; δODE, δcoarse) 

    NodesDict[:valleys]   = unique([to(γ) for γ in γ_to_valley]) # remove repeated valley if more than one contour goes there
    NodesDict[:entrances] = [to(γ) for γ in γ_to_entrance] 
    NodesDict[:poles] = unique([to(γ) for γ in γ_to_pole])

    plane_to_graph = Dict{ComplexF64,Int16}() # map complex plane points to graph vertices
    i = 0
    for key in keys(NodesDict)
        for z in NodesDict[key]
            i += 1
            plane_to_graph[z] = i
        end
    end

    ContourGraph = SimpleGraph(length(plane_to_graph))
    EdgesList = Dict{Tuple{Int,Int}, ComplexContour}()

    # Part 1: create edges between nodes inside the same ball.
    pts = [
        NodesDict[:exits]; NodesDict[:endpoint]; NodesDict[:entrances]; NodesDict[:statpoint]
    ]
    connect_inside_Ω!(pts, Ω, ContourGraph, plane_to_graph, EdgesList)

    # Part 2: create edges between stationary points with overlapping balls
    connect_stationarypoints!(Ω, ContourGraph, plane_to_graph, EdgesList)

    # Part 3: create edges between exit points and entrance points / valley points
    connect_ball_to_valleyyorentrance!(γ_to_entrance, ContourGraph, plane_to_graph, EdgesList)
    connect_ball_to_valleyyorentrance!(γ_to_valley, ContourGraph, plane_to_graph, EdgesList)
    connect_ball_to_valleyyorentrance!(γ_to_pole, ContourGraph, plane_to_graph, EdgesList)
    
    return ContourGraph, plane_to_graph, NodesDict, EdgesList
end

""" decide where to trace SD contours from"""

function _filter_exits(::AbstractPhaseFunction, Ω, exits)
    return exits # do nothing
end

function _filter_exits(::SquareRootPhaseFunction, Ω, exits)
    # this method is used avoid tracing SD contour at exit points inside Ω 
    # endpts = [a,b]
    if isreal.(exits) == [true,true] 
        return exits
        # if we moved along real line, check exits are outside non-osc region
        tracefrom = Vector{ComplexF64}()
        ex  = real.(exits) 
        δ = 1e-12 # by construction exits are exactly at the boundary, so we push them slightly away from Ω
        exδ = ex + sign.(ex .- Ω[1].c) * δ
        # if it remains in Ω after pushing them, then they are inside Ω and one should not trace from there
        @show !isinΩ(Ω,exδ[1]), !isinΩ(Ω,exδ[2])
        if !isinΩ(Ω,exδ[1]) push!(tracefrom, ex[1]) end
        if !isinΩ(Ω,exδ[2]) push!(tracefrom, ex[2]) end
        return tracefrom
    else
        return exits
    end
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
        inside_pts = [z for z in pts if abs(z - c) ≤ r + 10*eps()] 
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

function connect_ball_to_valleyyorentrance!(γvec :: Vector{ComplexContour}, CG ::SimpleGraph, CtoG::Dict, GtoContour::Dict)
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



function plot_ContourGraph(graph::SimpleGraph, Ω::Vector, z_to_G::Dict, NodesDict :: Dict)
    
    list, colors = _assign_colors_graph(z_to_G, NodesDict)

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
    # hidedecorations!(ax)
    ax.aspect = DataAspect()
    ax.backgroundcolor

    return fig

end

function plot_ContourGraph!(fig, ax, graph::SimpleGraph, Ω::Vector, z_to_G::Dict, NodesDict :: Dict)

    list, colors = _assign_colors_graph(z_to_G, NodesDict)
    _,_, p = graphplot!(ax, graph, node_color = colors, node_size = 20)
    mylayout(_) = list
    p.layout = mylayout

    # add non-oscillatory balls
    for Ball in Ω
        c,r = centre_and_radius(Ball)
        arc!(ax, Point2f(reim(c)), r, 0, 2π, 
            color = :gray)
    end

    hidespines!(ax)
    # hidedecorations!(ax)
    ax.aspect = DataAspect()
    ax.backgroundcolor

    return fig
end

function _assign_colors_graph(z_to_G::Dict, NodesDict::Dict)
    n = sum([length(NodesDict[k]) for k in keys(NodesDict)] )
    # Place graph nodes in the complex plane
    list = Vector(undef, n)

    # Add colors to nodes
    # Create a "Vector of colors" to pass an input to graph
    colors = Vector{Any}(undef, n)
    for z in NodesDict[:statpoint]
        i = z_to_G[z]
        colors[i] = colorant"red"
        list[i]   = reim(z)
    end
    for z in NodesDict[:endpoint]
        i = z_to_G[z]
        colors[i] = colorant"orange"
        list[i]   = reim(z)
    end
    for z in NodesDict[:exits]
        i = z_to_G[z]
        colors[i] = colorant"purple"
        list[i]   = reim(z)
    end
    R = _largest_point(NodesDict) + 1.0
    for z in NodesDict[:valleys]
        i = z_to_G[z]
        colors[i] = colorant"blue"
        list[i]   = reim(z) .* (R/abs(z))
    end
    for z in NodesDict[:entrances]
        i = z_to_G[z]
        colors[i] = colorant"green"
        list[i]   = reim(z)
    end
    for z in NodesDict[:poles]
        i = z_to_G[z]
        colors[i] = colorant"#5AC3E3"
        list[i]   = reim(z)
    end
    return list, colors
end

function _largest_point(NodesDict::Dict)
    R = 0.0
    for key in keys(NodesDict)
        if key != :valleys 
            if isempty(NodesDict[key]) continue end
            R = max(R, maximum(abs.(NodesDict[key])))
        end
    end
    return R
end