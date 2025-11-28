
function ContourGraph(G::AbstractPhaseFunction, a, b, Ω :: Vector{NonOscillatoryBall})
    # a,b are (finite) endpoints
    Pendp = [a,b]
    Pendp_trace = ComplexF64.([isinΩ(Ω,a) ? [] : a; isinΩ(Ω,b) ? [] : b])
    Pexit = exitpoints(G,Ω)
    Pstat = get_Pstat(Ω)

    global rstar = rvalley(G) # compute this only once

    SDpts = ComplexF64.([Pendp_trace; Pexit])
    # get entrance and valley points
    Pentr, Dict_entrance, Valleys, Dict_valleys  = tracing_contours(G, SDpts, Ω) 
    Valleys = unique(Valleys) # remove repeated valley if more than one contour goes there

    all_nodes = [z for z in [Pexit; [a,b]; Pentr; Valleys; Pstat]] 

    plane_to_graph = Dict{ComplexF64, Int16}() # map complex plane points to graph vertices
    for (i, z) in enumerate(all_nodes)
        plane_to_graph[z] = i
    end

    nvertices = length(plane_to_graph)
    ContourGraph = SimpleGraph(nvertices)

    # Part 1: create edges between nodes inside the same ball.
    pts = ComplexF64.([Pexit; Pendp; Pentr; Pstat])
    connect_inside_Ω!(pts, Ω, ContourGraph, plane_to_graph)

    # Part 2: create edges between stationary points with overlapping balls
    connect_stationarypoints!(Ω, ContourGraph, plane_to_graph)

    # Part 3: create edges between exit points and entrance points / valley points
    connect_ball_to_valleyyorentrance!(SDpts, Dict_entrance, ContourGraph, plane_to_graph)
    connect_ball_to_valleyyorentrance!(SDpts, Dict_valleys, ContourGraph, plane_to_graph)
    
    # @warn "SD contours going into the same valley are not connected"
    
    MetaDict = Dict{Symbol, Vector{ComplexF64}}()
    MetaDict[:valleys]   = Valleys
    MetaDict[:entrances] = Pentr
    MetaDict[:exits]     = Pexit
    MetaDict[:statpoint] = Pstat
    MetaDict[:endpoint]  = Pendp

    return ContourGraph, plane_to_graph, all_nodes, MetaDict
end

function connect_inside_Ω!(pts::Vector{ComplexF64}, Ω::Vector{NonOscillatoryBall}, CG ::SimpleGraph, dict::Dict)
    # create edges between nodes inside the same ball.
    # pts = [Pexit; Pendp; Pentr; Pstat]
    # CG is ContourGraph

    for Ball in Ω
        c, r = centre_and_radius(Ball)
        inside_pts = [z for z in pts if abs(z - c) ≤ r + 1e-14] 
        # added small perturbation here, sometimes the points at the boundary
        # are placed slightly outside the non-oscillatory region.

        for z1 in inside_pts
            for z2 in inside_pts
                if z1 != z2 add_edge!(CG, dict[z1], dict[z2]) end
            end
        end
    end

    return
end

function connect_stationarypoints!(Ω::Vector{NonOscillatoryBall}, CG ::SimpleGraph, dict::Dict)
    # create edges between stationary points if their respective balls overlap
    # CG is ContourGraph

    for Ball1 in Ω
        c1,r1 = centre_and_radius(Ball1)
        for Ball2 in Ω
            c2,r2 = centre_and_radius(Ball2)
            if Ball1 != Ball2 && (abs(c1-c2) < r1 +r2)
                add_edge!(CG, dict[c1], dict[c2])
            end
        end
    end
    return
end

function connect_ball_to_valleyyorentrance!(pts::Vector{ComplexF64}, ηdict::Dict, CG ::SimpleGraph, dict::Dict)
    # create edges between exit points - endpoints and entrance-valleys.
    # Dict contains η ∈ [Pexit, Pendp] and where it goes
    for η in pts
        if haskey(ηdict, η)
            goesto = ηdict[η]
            add_edge!(CG, dict[η], dict[goesto])
        end
    end
    return
end


function plot_ContourGraph(graph::SimpleGraph, Ω::Vector, nodes :: Vector{ComplexF64}, z_to_G::Dict, metadict :: Dict)
    
    # Place graph nodes in the complex plane
    list = Vector(undef, length(nodes))

    # Add colors to nodes
    # Create a "Vector of colors" to pass an input to graph
    colors = Vector{Any}(undef, length(nodes))
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