using PathFinder

ω    = 100000
f(z) = 1.0

# automated version
example1 = [0,0,0,0,0,0,0,0,0,1] # [3,5,6,2,9,5,1,4,1,3]

Phase = PolynomialPhaseFunction(example1) 
a,b = (-1,1) # specify (finite) endpoints 
val, figs = integrate(a,b,f,Phase,ω; 
                      plot_graph = true,
                      plot_sd = true)
figs[1]
figs[2]

# Airy phase debugging
x = -1
Phase = PolynomialPhaseFunction(-im*[0,-x,0,1/3])
a,b = (-π/3, π/3) # specify (infinite) endpoints 
val, figs = integrate(a,b,f,Phase,10.0; infcontour = [true,true],
                      plot_graph = true,
                      plot_sd = true)
figs[1]
figs[2]

# coalescence example debugging
P = 2
r = 1.0
coeffs = [0.0; -r^P; zeros(P-1); 1.0/(P+1)]  
Phase = PolynomialPhaseFunction(coeffs)
@show val,figs = integrate(-1,1,f,Phase,100.0; plot_graph = true, plot_sd = true)
@assert abs(-0.178593013067157 + 0.000000000000000im - val) < 1e-15
figs[1]
figs[2]