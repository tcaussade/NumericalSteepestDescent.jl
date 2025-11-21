
function ContourGraph(G::AbstractPhaseFunction, a, b, Ω :: Vector{NonOscillatoryBall})
    # a,b are (finite) endpoints
    Pendp = ComplexF64.([isinΩ(Ω,a) ? [] : a; isinΩ(Ω,b) ? [] : b])
    Pexit = exitpoints(G,Ω)
    Pstat = get_Pstat(Ω)


    global rstar = rvalley(G) # compute this only once

    SDpts = ComplexF64.([Pendp; Pexit])
    # get entrance and valley points
    Pentr, Dict_entrance, Valleys, Dict_valleys  = tracing_contours(G, SDpts, Ω) 

    all_nodes = [z for z in [Pexit; Pendp; Pentr; Valleys; Pstat]] 

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
    
    @warn "SD contours going into the same valley are not connected"

    return ContourGraph, plane_to_graph, all_nodes
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


function plot_ContourGraph(graph::SimpleGraph, nodes :: Vector{ComplexF64}, dict :: Dict)
    fig,ax,p = graphplot(graph)
    function mylayout(_)
        list = []
        for z in nodes
            push!(list, (real(z), imag(z)))
        end
        return list
    end

    p.layout = mylayout

    return fig

end