using PathFinder

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
acoefs = [0., 0, 1]
ps  = [0.]

ω = 20.0
a,b = -1im, 2



RatPhase = PathFinder.RationalPhaseFunction(acoefs, ps)
@show Ω = PathFinder.NonOscillatoryRegion(RatPhase, ω; Cball = 2π, δball = 1e-3,  Nrays = 16)
@show Pexit = PathFinder.exitpoints(RatPhase, Ω)

CG, CtoG, NodesDict, EdgesList = PathFinder.ContourGraph(RatPhase, a, b, Ω; δODE=0.1, δcoarse=0.01)

val, figs = integrate(a,b,z->1.0,RatPhase,ω; plot_graph = true, plot_sd = true)
figs[1]
figs[2]

using QuadGK

g(z) = PathFinder.evalphase(RatPhase, z)
refval, _ = quadgk(z-> cis(ω * g(z)), a,-1,im,b)
@show abs(refval - val) / abs(refval)


"""
    Questions and pendings:
    - non-oscillatory region might contain poles for low frequencies?
    - evaluation of residues (this is the hardest part!)
    - computing rstar for poles and valleys (currently it is an arbitrary value)
"""
