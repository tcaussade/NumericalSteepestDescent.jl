using PathFinder
using QuadGK

ω = 20.0
a,b = -2, 2

""" 
    Basic usage

acoefs = [c1,c2] is equivalent to an analytic part of c1 + c2*z
ps = [p] means there is a pole at p
pcoefs = [[d1,d2]] means that the pole a ps[1] takes the form d1/(z-p) + d2/(z-p)^2 
"""

# trying g(z) = z + 1/z
acoefs = [0,1]
ps     = [0.0]
pcoefs = [[1.]] #


RatPhase = RationalPhaseFunction(acoefs, ps, pcoefs)
# Ω = PathFinder.NonOscillatoryRegion(RatPhase, ω; Cball = 2π, δball = 1e-3,  Nrays = 16)
# Pexit = PathFinder.exitpoints(RatPhase, Ω)
# CG, CtoG, NodesDict, EdgesList = PathFinder.ContourGraph(RatPhase, a, b, Ω; δODE=0.1, δcoarse=0.01)
val, figs = PathFinder.integrate(a,b,z->1.0,RatPhase,ω; 
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]

Kp = 1
θ = π/(4Kp)
r = RatPhase.rstar_pole[1]
PathFinder.evaluate_noreturn_Gpole(r, θ, RatPhase; pole_idx = 1)

refval = quadgk(z -> cis(ω*(z+1/z)), a, im, b)[1]
@show abs(val-refval)/abs(refval)

"""
    Replicating Fig.3 in https://p-lpi.github.io

There, the phase is g(z) = 0.5*(x-μ)^2 + α/(1+x^2) with μ=0 and α=2
"""
ν = 200
α = 2
μ = 0

acoefs = [0.5*μ^2, -μ, 0.5]
ps = [im, -im]
pcoefs = 0.5im*α * [[-1], [1]]

PlasmaLensPhase = RationalPhaseFunction(acoefs, ps, pcoefs)
val, figs = PathFinder.integrate(π,0,z->1.0,PlasmaLensPhase,ν; 
                                infcontour = [true, true],
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]


idx = 1
Kp = 1
θ = π/(4Kp)
r = PlasmaLensPhase.rstar_pole[idx]
PathFinder.evaluate_noreturn_Gpole(r, θ, PlasmaLensPhase; pole_idx = idx)


""" 
    debugging catastrophe integral 
    See 36.2.6 in DLMF
"""
x,y,z = (2.0, 0.22, 0.0)
x,y,z = (5,5,0)
acoefs = [0,0,(x+z^2),0,2z,0,1]
ps     = [0.0]
pcoefs = [[0.0, y^2/12]] 
CatPhase = RationalPhaseFunction(acoefs, ps, pcoefs)

a,b = (-7π/12, π/12)
val, figs = PathFinder.integrate(a,b,z->1.0,CatPhase,100.0; infcontour = [true,true],
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]


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