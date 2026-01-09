using PathFinder

ω    = 10
f(z) = 1.0

# polynomial phase
example1 = [3,5,6,2,9,5,1,4,1,3]
Phase = PolynomialPhaseFunction(example1) 
z0,z1 = (-5,1) # specify (finite) endpoints 
val, figs = integrate(z0,z1,f,Phase,ω; 
                      plot_graph = true,
                      plot_sd = true)
figs[1]
figs[2]


# square-root phase g(z) = √(z^2+a^2) + bz
a,b = (0.00001, -0.7)
@show SqrtPhase = SquareRootPhaseFunction(a, b)
val, figs = integrate(z0, z1,f,SqrtPhase,ω; 
                      plot_graph = true,
                      plot_sd = true)
figs[1]
figs[2]
