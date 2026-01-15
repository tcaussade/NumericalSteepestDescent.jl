using PathFinder


# Basic usage example for arbitrary polynomial phase
ω    = 100 # frequency parameter
f(z) = 1.0 # Amplitude function
z0,z1 = (-1,1) # specify (finite) endpoints 

example1 = [3,5,6,2,9,5,1,4,1,3]
PolyPhase = PolynomialPhaseFunction(example1) 
val, figs = integrate(z0,z1,f,PolyPhase,ω; plot_graph = true, plot_sd = true)
figs[1]
figs[2]

# timings: MATLAB version takes around 30ms to evaluate with phase as above
@time integrate(z0,z1,f,PolyPhase,ω; plot_graph = false, plot_sd = false) # should take ~20ms 

# When the phase is linear the algorithm simplifies dramatically
LinPhase = LinearPhaseFunction()
val, _ = integrate(z0,z1, f, LinPhase, ω)

# We also can handle a square-root phase given by g(z) = √(z^2+a^2) + bz 
# This is a common integral to evaluate in HNA methods for high-frequency scattering problems
a,b = (1e-6, -0.0)
a,b = (0.35938136638046275, -0.1111111111111111)
SqrtPhase = SquareRootPhaseFunction(a, b)
@time val0, _ = integrate(0.0, 1.0,f,SqrtPhase,20; quadtype = :gaussian, N=20)

@time val1, figs = integrate(0.0, 1.0,f,SqrtPhase,20; quadtype = :adaptive, atol = 1e-8, plot_sd = true)
@show abs(val0 - val1)
