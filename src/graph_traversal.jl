
"""
    - For phases with no poles, it suffices to choose shortest path in the graph.

    - For rational phase, we might have to take a longer path to avoid resiudes.
    For this, we are using a Depth First Search (DFS) algorithm to find all contours connecting endpoints. 
"""

function get_deformation(::AbstractPhase, CG, start_point, end_point, CtoG, ::Dict, ::Dict, γ0)
    sd_edges = a_star(CG, CtoG[start_point], CtoG[end_point])
    return [(e.src, e.dst) for e in sd_edges]
end

# function quasi_sd_contour(::AbstractPhaseFunction, ::Dict, paths, ::Any)
#     idx = argmin(length.(paths))
#     return paths[idx]
# end


"""
    If G is a rational function,
    Find shortest path in the quasi SD contour that avoids poles.
"""

function get_deformation(G::RationalPhase, CG, start_point, end_point, CtoG, NodesDict::Dict, EdgesList::Dict, γ0)
    # extract nodes associated with valleys and poles
    vinf_nodes  = [CtoG[v] for v in NodesDict[:valleys]]
    vpole_nodes = [CtoG[v] for v in NodesDict[:poles]]
    vnodes = [vinf_nodes; vpole_nodes]

    all_SDedges = get_all_paths(CG, CtoG[start_point], CtoG[end_point], vnodes)
    sd_edges = quasi_sd_contour(G, EdgesList, all_SDedges, γ0)
    return sd_edges
end

function get_all_paths(CG::Graph, n1, n2, vnodes)

    paths_nodes = simple_paths(CG, n1, n2, vnodes)
    # paths_nodes = yen_k_shortest_paths(CG, n1, n2, weights(CG), 10).paths
    # @show length(paths_nodes)
    all_lists = []
    for path in paths_nodes
        edgelist = [(path[i], path[i+1]) for i in eachindex(path[1:end-1])]
        push!(all_lists, edgelist)
    end
    return all_lists
end

function quasi_sd_contour(G::RationalPhase, EdgesList::Dict, paths, γ0)
    # γ0 is the collection of nodes for starting integration contour.
    sort!(paths, by = length) # sort by length and begin with shortest
    for path in paths
        γ = Vector{ComplexContour}()
        for e in path # get a quasi-sd contour
            haskey(EdgesList,e) ? push!(γ, EdgesList[e]) : push!(γ, EdgesList[reverse(e)])
        end
        nγ = zero(ComplexF64)
        for zp in poles(G) # check if γ - γ0 crossed poles
            γfull = winding_contour(γ, γ0)
            nγ += winding_number(zp,γfull)
        end
        # TOLERANCE HERE SHOULD BE MORE STRICT?
        if abs(nγ)<0.1 return path end
    end
    @warn "There are no residue-free paths"
    return
end



"""
    Custom implementation of a graph traversal algorithm to find paths between two nodes in a graph,
     with the additional constraint that the path should not cross any poles of the integrand. 
"""

function simple_paths(G::Graph, n1, n2, valley_nodes :: Vector{Int16})
    # This function returns all simple paths between n1 and n2 in the graph G
    # We can use a Depth First Search (DFS) algorithm to find all paths
    # Warning: this algorithm has factorial complexity with respect to the number of nodes, so it is not suitable for large graphs.
    visited = falses(nv(G))
    path = Int[]
    all_paths = Vector{Vector{Int}}()
    
    function dfs(v)
        push!(path, v)
        # we allow to visit nodes associated valley_nodes multiple times, but not other nodes
        if v in valley_nodes
            visited[v] = false
        else
            visited[v] = true
        end

        if v == n2
            push!(all_paths, copy(path))
        else
            for w in neighbors(G, v)
                if !visited[w]
                    dfs(w)
                end
            end
        end
        
        pop!(path) # reset visited status for backtracking
        visited[v] = false
    end
    
    dfs(n1)
    return all_paths
end

# # remove redundant edges in paths that are traversed in both directions
# function filter_path(path)
#     fil

# end