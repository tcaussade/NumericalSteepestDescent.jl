using Graphs
using GraphMakie

# Goal: Create a graph with 5 vertices and random edges
# we should be able to plot it later, with nodes placesd in the compelx plane

n = 4
z = rand(ComplexF64, n)
z = [1+im, 2+ 1.2im, -3, -1im]
z_to_graph = Dict()
for (i,z) in enumerate(z)
    z_to_graph[z] = i
end

G = SimpleGraph(n)
for z1 in z
    for z2 in z
        if z1 != z2
            println("Adding edge between $z1 and $z2")
            # add edge to graph
            add_edge!(G, z_to_graph[z1], z_to_graph[z2])
        end
    end
end

fig,ax,p = graphplot(G)

p.layout = mylayout
fig


function mylayout(_)
    list = []
    for z in z
        push!(list, (real(z), imag(z)))
    end
    return list
end