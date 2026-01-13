using PathFinder

ω    = 10
f(z) = 1.0
z0,z1 = (-1,1) # specify (finite) endpoints 

# polynomial phase
example1 = [3,5,6,2,9,5,1,4,1,3]
PolyPhase = PolynomialPhaseFunction(example1) 
val, figs = integrate(z0,z1,f,PolyPhase,ω; plot_graph = true, plot_sd = false)
figs[1]
figs[2]

# linear phase 
LinPhase = LinearPhaseFunction()
val, figs = integrate(z0,z1, f, LinPhase, ω; plot_sd = true, plot_graph = true)
figs[2]

# square-root phase g(z) = √(z^2+a^2) + bz
a,b = (1., 0.)
SqrtPhase = SquareRootPhaseFunction(a, b)
val, figs = integrate(0.0, z1,f,SqrtPhase,ω; plot_sd = true, plot_graph = true)
figs[2]

# do tests! especially for b = ±1

x = 0.0
y = 0.0
G = PolynomialPhaseFunction([0.0, 0.0, 0.0, 1.0])
I_GHH,figs = integrate(π, 0.0, z -> 1.0, G, 1.0; N=20, infcontour = [true true], plot_sd = true, plot_graph = true)
figs[1]