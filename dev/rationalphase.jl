using PathFinder

# struct RationalPhaseFunction
#     num :: Polynomial # numerator
#     den :: Polynomial # denominator
#     ξ   :: Vector # stationary points
#     p   :: Vector # poles
#     function RationalPhaseFunction(num_coefs ::Vector,pole_vals::Vector)
#         num = Polynomial(num_coefs)
#         den = fromroots(pole_vals) # Polynomial(den_coefs)
#         dnum = derivative(num)*den - num*derivative(den)
#         ξ = roots(dnum)  
#         p = pole_vals
#         new(num, den, ξ ,p)
#     end
# end

# poles(R::RationalPhaseFunction) = R.p
# phasefunction(R :: RationalPhaseFunction) = R.num // R.den

# # z + 1/z = (z^2 + 1) / z
# num,ps = [1,0,1],[0]
# num,ps = [0,1],[im, -im]
# num,ps = [1], [1,1,-2]
# RatPhase = RationalPhaseFunction(num, ps)

# function evalResidue(R::RationalPhaseFunction, z0)
#     @assert z0 in poles(R)
#     @show n = count(==(z0), poles(R)) # count multiplicity
#     lim = derivative(R.num // fromroots(setdiff(poles(R), z0)), n-1)
#     return lim(z0) * 1/factorial(n-1)
# end

# evalResidue(RatPhase, 1)

""" trying g(z) = z + 1/z = (z^2+1)/z """
num = [1, 0, 1]
ps  = [0]

ω = 1.0
a,b = 1,2

RatPhase = PathFinder.RationalPhaseFunction(num, ps)
Ω = PathFinder.NonOscillatoryRegion(RatPhase, ω; Cball = 2π, δball = 1e-3,  Nrays = 16)
CG, CtoG, NodesDict, EdgesList = PathFinder.ContourGraph(RatPhase, a, b, Ω; δODE=0.1, δcoarse=0.01)