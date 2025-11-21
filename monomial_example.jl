using NumericalStationaryPhase

ω     = 40.0
Cball = 2π
J     = 2
f(z) = 1.0

Monomial = MonomialPhaseFunction(J)

# Quasi-SD contour deformation
Ω = NonOscillatoryRegion(Monomial, Cball, ω)
γ = QuasiSDcontour(Monomial, Ω[1], 1.0)

# plot contour deformation
plot_quasiSDdeformation(Monomial, γ, Ω)

# evaluate the integral
using FastGaussQuadrature
x1,w1 = gausslegendre(20);
x2,w2 = gausslaguerre(20);

I = 0.0 + 0im

c1 = γ[1] # finite contour
I += eval_finite(f, Monomial, c1, ω, x1, w1)

c2 = γ[2] # infinite SD contour at boundary of NonOscRegion
I += eval_infiniteSDpath(f, Monomial, c2, ω, x2, w2)
    
c3 = γ[3] # infinite SD contour at endpoint
I += -eval_infiniteSDpath(f, Monomial, c3, ω, x2, w2)
   

N = 20
In = integrate_nsp(f, Monomial, γ, ω, N) ≈ I # should be true
