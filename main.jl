using PathFinder

# Basic usage example for arbitrary polynomial phase
ω    = 10 # frequency parameter
f(z) = 1.0 # Amplitude function
z0,z1 = (-1,1) # specify (finite) endpoints 

example1 = [3,5,6,2,9,5,1,4,1,3]
example1 = ([0,+0.05,0,1])
example1 = [0,0,1]
PolyPhase = PolynomialPhaseFunction(example1) 
val, figs = integrate(-1,1,f,PolyPhase,ω; plot_graph = true, plot_sd = true)
figs[1]
figs[2]

_,fig = integrate(π, 0.0, x -> 1.0, PolynomialPhaseFunction([0, -1, +0.0, 0, 1]), 1.0, infcontour=[true,true],
plot_sd = true, Nrays = 5) # should take ~0.3ms with @benchmark 

bm = @benchmark val = integrate(-1,1,f,PolyPhase,ω)

# When the phase is linear the algorithm simplifies dramatically
LinPhase = LinearPhaseFunction()

""" quick test for polynomial """
r = PolyPhase.rstar_valley
J = length(example1)-1
PathFinder.evaluate_noreturn_Ginf(r,π/(4J),PolyPhase)

# When the phase is linear the algorithm simplifies dramatically
LinPhase = LinearPhaseFunction()
val, _ = integrate(z0,z1, f, LinPhase, ω)

# We also can handle a square-root phase given by g(z) = √(z^2+a^2) + bz 
# This is a common integral to evaluate in HNA methods for high-frequency scattering problems
a,b = (0.1, 0.0)
ω = 50
SqrtPhase = SquareRootPhaseFunction(a, b)
val0, figs = integrate(0.0, 1.0,f,SqrtPhase,ω; quadtype = :gaussian, N=10, 
                plot_sd=true, plot_graph=true)

figs[1]
figs[2]

# val1, figs = integrate(0.0, 1.0,f,SqrtPhase,20; quadtype = :adaptive, atol = 1e-8, plot_sd = true)
# @show abs(val0 - val1)
