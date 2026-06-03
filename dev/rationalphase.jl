using PathFinder
using QuadGK

ω = 20

""" 
    Basic usage

acoefs = [c1,c2] is equivalent to an analytic part of c1 + c2*z
ps = [p] means there is a pole at p
pcoefs = [[d1,d2]] means that the pole a ps[1] takes the form d1/(z-p) + d2/(z-p)^2 
"""

# trying g(z) = z + 1/z
acoefs = [0,1]
ps     = [0.0]
pcoefs = [[1]] #


RatPhase = RationalPhase(acoefs, ps, pcoefs)
# Ω = PathFinder.NonOscillatoryRegion(RatPhase, ω; Cball = 2π, δball = 1e-3,  Nrays = 16)
# Pexit = PathFinder.exitpoints(RatPhase, Ω)
# CG, CtoG, NodesDict, EdgesList = PathFinder.ContourGraph(RatPhase, a, b, Ω; δODE=0.1, δcoarse=0.01)
val, figs = PathFinder.nsd([-2,-im,2],z->1.0,RatPhase,ω; 
                                infcontour = [false, false],
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]

nsd([-2,-1im,2],z->1.0,RatPhase,ω)

v= NodesDict[:valleys]
vnodes = [CtoG[v] for v in NodesDict[:valleys]]

p1 = yen_k_shortest_paths(CG, 1, 2, weights(CG), 10).paths
p2 = all_simple_paths(CG, Int16(1), Int16(2)) |> collect
p3 = simple_paths(CG, Int16(1), Int16(2),vnodes)

# bm = @benchmark integrate(-2,2,z->1.0,RatPhase,ω)

Kp = 1
θ = π/(4Kp)
r = RatPhase.rstar_pole[1]
PathFinder.evaluate_noreturn_Gpole(r, θ, RatPhase; pole_idx = 1)

r = RatPhase.rstar_valley
J = 2
θ = π/(4J)
PathFinder.evaluate_noreturn_Ginf(r,θ, RatPhase)

refval = quadgk(z -> cis(ω*(z+1/z)), a, im, b)[1]
@show abs(val-refval)/abs(refval)

"""
    Replicating Fig.3 in https://p-lpi.github.io

There, the phase is g(z) = 0.5*(x-μ)^2 + α/(1+x^2) with μ=0 and α=2
"""
ν = 20
μ, α = (-0.4,1.0344893103448276)
μ,α = (0,5)

acoefs = [0.5*μ^2, -μ, 0.5]
ps = [im, -im]
pcoefs = 0.5im*α * [[-1], [1]]

PlasmaLensPhase = RationalPhase(acoefs, ps, pcoefs)
val, figs = PathFinder.nsd([π,0],z->1.0,PlasmaLensPhase,ν; 
                                infcontour = [true, true],
                                plot_graph = true, plot_sd = true)
figs[2]

# r = PlasmaLensPhase.rstar_valley
# J = 2
# θ = π/(4J)
# PathFinder.evaluate_noreturn_Ginf(r,θ, PlasmaLensPhase)

# idx = 1
# Kp = 1
# θ = π/(Kp) * 0.2500001
# r = PlasmaLensPhase.rstar_pole[idx]
# PathFinder.evaluate_noreturn_Gpole(r, θ, PlasmaLensPhase; pole_idx = idx)


PathFinder.winding_number(-0.9999im, [-1.0-im, 1.0-im, 1.0+im, -1.0+im, -1.0-im]) # should be 1




""" 
    debugging catastrophe integral 
    See 36.2.6 in DLMF
"""
x,y,z = (2.0, 0.22, 0.0)
x,y,z = (0,5,0)
x,y,z = (3.673469387755102,-1.7142857142857142,0.0)
x,y,z = (3.6,-1,0.0)
acoefs = [0,0,(x+z^2),0,2z,0,1]
ps     = [0.0]
pcoefs = [[0.0, y^2/12]] 
CatPhase = RationalPhase(acoefs, ps, pcoefs)

a,b = (-7π/12, π/12)
val, figs = PathFinder.nsd([a,b],z->1.0,CatPhase,1; infcontour = [true,true],
                                plot_graph = true, plot_sd = true)
figs[1]
figs[2]

@profview nsd(a,b,z->1.0,CatPhase,1; infcontour = [true,true])

""" arbitrary """
d = 0.5
Dipole = RationalPhase([3,1,3,4,5,1,2], [d,-d], [[0.5,-0.25im],[0.5,0.5,0.5]])
v,fig = nsd([-im,im,1],z->1.0,Dipole,10.0; plot_graph = true, plot_sd = true)
fig[2]

""" Case when poles might lie close to each other """

acoefs = [0,0,1]
ps     = [0.0, 0.01, 0.5, 1.0,1.1]
pcoefs = [[1], [-im], [1], [1], [1]]
ClosePhase = RationalPhase(acoefs, ps, pcoefs)
val, figs = PathFinder.nsd(-1,2,z->1.0,ClosePhase,10.0;
                                plot_graph = true, plot_sd = true)
figs[2]

""" Arbitrary Rational phase ? """

acoefs = [3,1,4]
ps     = [0, im]
pcoefs = [[1,2,3], [im]]# [[1,6,1,8], [4,9], [4]]
RatPhase = RationalPhase(acoefs, ps, pcoefs)

val,plt = PathFinder.nsd(-1,1,z->1.0,RatPhase,10.0; plot_sd = true)

fig = PathFinder.Figure()
ax  = PathFinder.Axis(fig[1, 1], title = "", aspect = PathFinder.DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -2:1:2, yticks = -2:1:2)
PathFinder.quasiSDdeformation!(fig, ax, -1,1, RatPhase, 1;
                                                    umax = 80,
                                                    color_lim = 300,
                                                    resolution = 200,
                                                    set = 8)

"""
    Questions and pendings:
    - non-oscillatory region might contain poles for low frequencies?
    - evaluation of residues (this is the hardest part!)
    - when stationary point lie close to poles, they might stay inside the rstar_pole ball, which leads to missing edges in graph
    - sometimes Pexit fails to find all exit points - can we devise an algorithm that always finds the solutions?
    - 
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

