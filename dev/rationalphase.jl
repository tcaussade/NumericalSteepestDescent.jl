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
acoefs = im*[0,0,-2,0,0,0,1] # [0,1]
ps     = [0.0]
pcoefs = [[0.0, im*1/3]] # Should be a Vector of vectors to allow different pole orders at different points
pcoefs = 

ω = 1.0
a,b = -1, 2

RatPhase = PathFinder.RationalPhaseFunction(acoefs, ps, pcoefs)
@show RatPhase.vpole/π
Ω = PathFinder.NonOscillatoryRegion(RatPhase, ω; Cball = 2π, δball = 1e-3,  Nrays = 16)
Pexit = PathFinder.exitpoints(RatPhase, Ω)
CG, CtoG, NodesDict, EdgesList = PathFinder.ContourGraph(RatPhase, a, b, Ω; δODE=0.1, δcoarse=0.01)
val, figs = PathFinder.integrate(a,b,z->1.0,RatPhase,ω; 
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]




a,b = (-7π/12, π/12) .-π/24
val, figs = PathFinder.integrate(a,b,z->1.0,RatPhase,ω; infcontour = [true,true],
                                plot_graph = false, plot_sd = true)
figs[2]

using QuadGK

g(z) = PathFinder.evalphase(RatPhase, z)
# refval, _ = quadgk(z-> cis(ω * g(z)), a,-1,im,b)
refval,_ = quadgk(z->cis(ω*g(z)), a,b)
@show abs(refval - val) / abs(refval)


"""
    Questions and pendings:
    - non-oscillatory region might contain poles for low frequencies?
    - evaluation of residues (this is the hardest part!)
    - computing rstar for poles and valleys (currently it is an arbitrary value)
    - the root finding methods struggle to converge!!
"""

"""
    Ideas for residues:
        Maybe use quadgk to evaluate them on small discs around poles?
"""

"""Order to do stuff
        1. compute rstar correctly
        2. think about ball radius, should I change something?
        3. residues?
"""

acoefs = [0,0,1]
poles  = [0.0]
poles_coefs = [[0.0 1.0];]
function def_ratphase(acoefs, poles, poles_coefs)
    @show apart = Polynomial(acoefs)
    id = Polynomial(1.0)
    spart = Polynomial(0.0)
    for (i,zp) in enumerate(poles)
        for (k,coef) in enumerate(poles_coefs[i,:])
            pvec = zp * ones(k) 
            spart += coef * id // fromroots(pvec) 
        end
    end
    @show spart
    return lowest_terms(apart + spart)
end

def_ratphase(acoefs, poles, poles_coefs)