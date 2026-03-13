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

# remove redundant edges in paths that are traversed in both directions
function filter_path(path)
    fil

end