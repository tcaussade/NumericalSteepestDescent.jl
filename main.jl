using PathFinder
# using Polynomials

ω    = 1.0
f(z) = 1.0


x = -1.1
Phase = PolynomialPhaseFunction(-im*[0,-x,0,1/3])

# Phase = PolynomialPhaseFunction([0,0,0,0,1])

example1 = [3,5,6,2,9,5,1,4,1,3]
Phase = PolynomialPhaseFunction(example1) 

# Phase = PolynomialPhaseFunction([0,0,1])

# automated version
a,b = (-1-im,2+2im) # specify (finite) endpoints
a,b = (-π/3, π/3) # specify (infinite) endpoints 
val, figs = integrate(a,b,f,Phase,ω; 
                    #   infcontour = [true,true],
                      plot_graph = true,
                      plot_sd = true)

figs[1]
figs[2]

P = 6
r = 1.0
coeffs = [0.0; -r^P; zeros(P-1); 1.0/(P+1)]  
Phase = PolynomialPhaseFunction(coeffs)
_,figs = integrate(-1,1,f,Phase,100.0; plot_graph = true, plot_sd = false)
figs[2]
figs[1]