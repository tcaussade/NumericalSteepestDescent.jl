using PathFinder
# using Polynomials

ω     = 40.0
Cball = 2π
f(z) = 1.0

example1 = [3,5,6,2,9,5,1,4,1,3]

acubic = -1
Phase = PolynomialPhaseFunction([0,-3*acubic,0,1])
# Phase = PolynomialPhaseFunction(example1) 

# Quasi-SD contour deformation
a,b = (0,2) # specify (finite) endpoints
Ω = NonOscillatoryRegion(Phase, Cball, ω)

Pexit = PathFinder.exitpoints(Phase,Ω)
Pstat = PathFinder.get_Pstat(Ω)

graph, dict, nodes, metadict = PathFinder.ContourGraph(Phase, a,b, Ω)
PathFinder.plot_ContourGraph(graph,Ω, nodes, dict, metadict)

# get shortest path
path = a_star(graph, dict[a], dict[b])
path[1].src, path[1].dst # source and destination of edge 1

# automated version
γ = QuasiSDcontour(Phase, a,b, Ω)


# plot contour deformation
plot_quasiSDdeformation(Monomial, γ, Ω)



