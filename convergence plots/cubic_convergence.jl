using NumericalSteepestDescent
using CairoMakie
# using WGLMakie
# using Polynomials

# choose roots
# dP = fromroots([-1.1,0.3im, 1+im, 1.2-im/5])
# P = Polynomials.integrate(dP)

freqs  = ([25, 50,100, 200, 400, 800])
# freqs  = ([50, 200, 500])
# freqs = [20,40,80,160,320,1280]
# freqs = [20,40,80,100]
nquads = range(1,25) 

a = -1/3 * 1
CubicPhase = PolynomialPhaseFunction(-im*[0,-3a,0,1])
# CubicPhase = PolynomialPhaseFunction([3,5,6,2,9,5,1,4,1,3])
# CubicPhase = PolynomialPhaseFunction(P.coeffs)
# CubicPhase = PolynomialPhaseFunction([0,1,2,3,4,5,6] * 0.1)
# x,y = 1,2
# CubicPhase = PolynomialPhaseFunction([0, x, y, 0, 1])

# P = 6
# r = 0.1
# CubicPhase = PolynomialPhaseFunction([0.0; -r^P; zeros(P-1); 1.0/(P+1)])

# CubicPhase = PolynomialPhaseFunction([0, -8+24im, -(10+44im)/3, -(4-28im)/9, (1+im), -4/15, 1/18])

val,fig = integrate([-1,1],z->1.0,CubicPhase,100; plot_sd = true)
@show val
fig[1]

f(z) = 1 #/(z-1-0.1im)
e = zeros(length(freqs), length(nquads))
for (i,ω) in enumerate(freqs)
    ref = PathFinder.integrate([-π/3,π/3],f,CubicPhase,ω; N = 50, infcontour = [true, true])
    for (n,N) in enumerate(nquads)
        val = PathFinder.integrate([-π/3,π/3],f,CubicPhase,ω; N, infcontour = [true true], Cball = 8π/8)
        e[i,n] = abs(ref.-val) / abs(ref)
    end
end

fig = Figure()
ax = Axis(fig[1,1], 
            xlabel = "number of quadrature points per contour (sqrt scale)",
            ylabel = "absolute error (log scale)",
            yscale = log10, xscale = sqrt)

for (i,ω) in enumerate(freqs)
    scatterlines!(nquads, e[i,:], label = "ω = $ω")
end
axislegend(ax)
limits!(1,26,1e-16,1e-2)
fig

# uniformity in parameter space

nxPts = 500
xArray = range(-8, 5, length=nxPts)

Err = zeros(nxPts)
fig = PathFinder.Figure()
ax = PathFinder.Axis(fig[1,1], 
            xlabel = "x",
            ylabel = "absolute error (log scale)",
            yscale = log10)

for nquad in [5,10,15,20,30]
    for (i,x) in enumerate(xArray)
        airyJulia = airyai(x)
        G = PolynomialPhaseFunction(-1im * [0, -x, 0, 1/3])
        integral = PathFinder.integrate(-π/3, π/3, x -> 1.0, G, 1.0; N=nquad, infcontour=[true,true])
        airyPathFinder = (1/(2im*π)) * integral
        Err[i] = abs(airyJulia - airyPathFinder) #/ abs(airyJulia)
    end
    PathFinder.scatterlines!(xArray, Err, label = "N = $nquad")
end
PathFinder.axislegend(ax, position = :rc)
PathFinder.limits!(-8.1,5.1,1e-18,1e-0)
fig

