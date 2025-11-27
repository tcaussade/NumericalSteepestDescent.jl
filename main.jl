using NumericalStationaryPhase
# using Polynomials

ω     = 40.0
Cball = 2π
f(z) = 1.0

# example1 = [3,5,6,2,9,5,1,4,1,3]

acubic = -1
Phase = PolynomialPhaseFunction([0,-3*acubic,0,1])
# Phase = PolynomialPhaseFunction(example1) 

NumericalStationaryPhase.goes_to_valley(Phase, )

# Quasi-SD contour deformation


a,b = (0,1) # specify (finite) endpoints
Ω = NonOscillatoryRegion(Phase, Cball, ω)

Pexit = NumericalStationaryPhase.exitpoints(Phase,Ω)

graph, dict, nodes, metadict = NumericalStationaryPhase.ContourGraph(Phase, a,b, Ω)
NumericalStationaryPhase.plot_ContourGraph(graph, nodes, dict, metadict)

# get shortest path
sd = dijkstra_shortest_paths(graph, dict[a])
sd.dists

# automated version
γ = QuasiSDcontour(Phase, a,b, Ω)


# plot contour deformation
plot_quasiSDdeformation(Monomial, γ, Ω)



