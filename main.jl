using PathFinder
# using Polynomials

ω     = 4.0
Cball = 2π
f(z) = 1.0

example1 = [3,5,6,2,9,5,1,4,1,3]

acubic = +0.2
Phase = PolynomialPhaseFunction([0,-3*acubic,0,1])
Phase = PolynomialPhaseFunction(example1) 
# Phase = PolynomialPhaseFunction([0,0,1])

# automated version
a,b = (-1,1) # specify (finite) endpoints
val = PathFinder.integrate(a,b,f,Phase,ω)


# Quasi-SD contour deformation

Ω = NonOscillatoryRegion(Phase, Cball, ω)

Pexit = PathFinder.exitpoints(Phase,Ω)
Pstat = PathFinder.get_Pstat(Ω)

graph, dict, metadict, edgeslist = PathFinder.ContourGraph(Phase, a,b, Ω)
PathFinder.plot_ContourGraph(graph,Ω,dict, metadict)
# there is a bug when endpoiints coincide with stationary points

# get shortest path
path = a_star(graph, dict[a], dict[b])
path[1].src, path[1].dst # source and destination of edge 1


# plot contour deformation
plot_quasiSDdeformation(Monomial, γ, Ω)



