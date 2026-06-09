using NumericalSteepestDescent
using QuadGK

a = -1
Cubic = PolynomialPhase([0,-3a,0,1])

ω = 1
L = 1e-12
i1 = nsd([0,L],x -> 1,Cubic,ω; N = 50)
i2 = quadgk(z -> cis(ω*Cubic.p(z)), 0, L)[1]
@show abs(i1-i2)/abs(i2)


####
using FastGaussQuadrature

x,w = gausslaguerre(20)
s = 100+im
f(x) = 1/(x-s)
i1 = sum(w .* f.(x))
i2 = quadgk(x -> f(x)*exp(-x),0,Inf)[1]
@show abs(i1-i2) #/abs(i2)


Mon = PolynomialPhase([0,0,1])
ω = 40
i1,fig = nsd([cis(-π/4),1],x -> 1,Mon,ω; 
                N = 50, infcontour = [false, false], 
                plot_sd = true)
i1
limits!(-2,2,-2,2)
fig[1]
