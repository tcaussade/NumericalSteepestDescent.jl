using PathFinder
# using Polynomials

ω     = 2.0
Cball = 2π
f(z) = 1.0

example1 = [3,5,6,2,9,5,1,4,1,3]

acubic = +0.2
Phase = PolynomialPhaseFunction([0,-3*acubic,0,1])
Phase = PolynomialPhaseFunction(example1) 
# Phase = PolynomialPhaseFunction([0,0,1])

# automated version
a,b = (-1,1) # specify (finite) endpoints
val, figs = integrate(a,b,f,Phase,ω; 
                      plot_graph = false,
                      plot_sd = true)
figs[2]

# Quasi-SD contour deformation

# Ω = NonOscillatoryRegion(Phase, Cball, ω)

# Pexit = PathFinder.exitpoints(Phase,Ω)
# Pstat = PathFinder.get_Pstat(Ω)

# graph, dict, metadict, edgeslist = PathFinder.ContourGraph(Phase, a,b, Ω)
# PathFinder.plot_ContourGraph(graph,Ω,dict, metadict)
# # there is a bug when endpoiints coincide with stationary points


