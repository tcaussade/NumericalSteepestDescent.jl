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
acoefs = [0,1]
ps     = [0.0]
pcoefs = [[0.0, 1.]] # Should be a Vector of vectors to allow different pole orders at different points

ω = 1.0
a,b = -2, 2

RatPhase = RationalPhaseFunction(acoefs, ps, pcoefs)
Ω = PathFinder.NonOscillatoryRegion(RatPhase, ω; Cball = 2π, δball = 1e-3,  Nrays = 16)
Pexit = PathFinder.exitpoints(RatPhase, Ω)
CG, CtoG, NodesDict, EdgesList = PathFinder.ContourGraph(RatPhase, a, b, Ω; δODE=0.1, δcoarse=0.01)
val, figs = PathFinder.integrate(a,b,z->1.0,RatPhase,ω; 
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]

""" trying catastrophe integral """
x,y,z = (2.0, 0.2222222, 0.0)
acoefs = im*[0,0,2(x+z^2),0,2z,0,1]
ps     = [0.0]
pcoefs = [[0.0, im*y^2/12]] 
RatPhase = RationalPhaseFunction(acoefs, ps, pcoefs)

a,b = (-7π/12, π/12) .-π/24
val, figs = PathFinder.integrate(a,b,z->1.0,RatPhase,5.0; infcontour = [true,true],
                                plot_graph = true, plot_sd = true)
figs[1]
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
    - when stationary point lie close to poles, they might stay inside the rstar_pole ball, which leads to missing edges in graph
    - sometimes Pexit fails to find all exit points.
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