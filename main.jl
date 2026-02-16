using PathFinder

# Basic usage example for arbitrary polynomial phase
ω    = 50 # frequency parameter
f(z) = 1.0 # Amplitude function
z0,z1 = (-1,1) # specify (finite) endpoints 

example1 = [3,5,6,2,9,5,1,4,1,3]
example1 = ([0,+0.05,0,1])
PolyPhase = PolynomialPhaseFunction(example1) 
val, figs = integrate(z0,z1,f,PolyPhase,ω; plot_graph = true, plot_sd = true)
figs[1]
figs[2]

using BenchmarkTools
# # timings: MATLAB version takes around 30ms to evaluate with phase as above
# @time integrate(z0,z1,f,PolyPhase,ω; plot_graph = false, plot_sd = false) # should take ~2ms with @benchmark 
# # MATLAB takes around 5ms per evaluation for cuspoid integral
@benchmark integrate(π, 0.0, x -> 1.0, PolynomialPhaseFunction([0, -0.01, +0.2, 0, 1]), 1.0, infcontour=[true,true]) # should take ~0.3ms with @benchmark 

# When the phase is linear the algorithm simplifies dramatically
LinPhase = LinearPhaseFunction()
val, _ = integrate(z0,z1, f, LinPhase, ω)

# We also can handle a square-root phase given by g(z) = √(z^2+a^2) + bz 
# This is a common integral to evaluate in HNA methods for high-frequency scattering problems
a,b = (1.0, -1/sqrt(2))
a,b = (1.0, -0.8947368421052632)
ω = 50
@show SqrtPhase = SquareRootPhaseFunction(a, b)
val0, figs = integrate(0.0, 1.0,f,SqrtPhase,ω; quadtype = :gaussian, N=10, 
                plot_sd=true, plot_graph=true)

figs[1]
figs[2]

# val1, figs = integrate(0.0, 1.0,f,SqrtPhase,20; quadtype = :adaptive, atol = 1e-8, plot_sd = true)
# @show abs(val0 - val1)
