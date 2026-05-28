using NumericalSteepestDescent
using QuadGK
# using WGLMakie
using SpecialFunctions

quadtype = :gaussian
quads = [10,15,30, 100]
# ωvals = [1,5,10,20]
ω = 20

singular = false

avals = 10 .^ range(-8, 0, length=80)
bfix  = 0.0
relErrA = zeros(length(avals), length(quads))

bvals = range(-1.0,1.0, length = 80)[2:end-1]
afix  = 1.0
relErrB = zeros(length(bvals), length(quads))

# for ω in ωvals
    
    for (j,a) in enumerate(avals)

        fA(z) = singular ? hankelh1(0, ω*sqrt(z^2 + a^2)) * cis(-ω * sqrt(z^2 + a^2)) : 1.0
        ζA(z) = fA(z) * cis(ω * (sqrt(z^2+a^2) + bfix*z))
        refA = quadgk(ζA, 0.0, 1.0, atol = 1e-14)[1] # brute-force
        G = SquareRootPhaseFunction(a,bfix)
        @show (a, bfix)
        for (n,N) in enumerate(quads)
            intA = integrate(0.0, 1.0, fA, G, ω; N)  
            relErrA[j,n] = abs(intA-refA) #/ abs(refA)      
        end   
    end
    
    
    for (j,b) in enumerate(bvals)
        fB(z) = singular ? hankelh1(0, ω*sqrt(z^2 + afix^2)) * cis(-ω * sqrt(z^2 + afix^2)) : 1.0
        ζB(z) = fB(z) * cis(ω * (sqrt(z^2+afix^2) + b*z))
        refB = quadgk(ζB, 0.0, 1.0, atol = 1e-14)[1] # brute-force
        G = SquareRootPhaseFunction(afix,b)
        @show (afix, b)
        for (n,N) in enumerate(quads)
            intB = integrate(0.0, 1.0, fB, G, ω; N)[1]   
            relErrB[j,n] = abs(intB-refB)#/ abs(refB)      
        end   
    end
          
# end

fig = PathFinder.Figure(size = (800,400))
axA = PathFinder.Axis(fig[1, 1], title = "Absolute error for a ∈ (0,1) with b = $bfix",  
                xlabel = "a (log scale)", ylabel = "Abs. error", limits = (nothing, nothing,1e-16,1e-0), yscale = log10, xscale = log10)
axB = PathFinder.Axis(fig[1, 2], title = "Absolute error for b ∈ [-1,1] with a = $afix",  
                xlabel = "b", ylabel = "Abs. error", limits = (-1,1,1e-16,1e-0), yscale = log10)
for (i,n) in enumerate(quads)
    PathFinder.lines!(axA, avals, relErrA[:,i], label = "n = $n")
    PathFinder.lines!(axB, bvals, relErrB[:,i], label = "n = $n")
end
PathFinder.axislegend(axA)
PathFinder.axislegend(axB)
fig

# Numerical convergence

avalsN = [0.1, 1]
nquads = 1:30
f(z) = 1 #/(z-1-0.1im)
ω = 50

e = zeros(length(avalsN), length(nquads))
for (i,a) in enumerate(avalsN)
    Phase = SquareRootPhaseFunction(a,0)    
    ref = PathFinder.integrate(0,1,f,Phase,ω; N = 100, Cball = 3pi) # Cball = 8π/8)
    for (n,N) in enumerate(nquads)
        val = PathFinder.integrate(0,1,f,Phase,ω; N, Cball = 3π)
        e[i,n] = abs(ref.-val) #/ abs(ref)
    end
end

fig = PathFinder.Figure()
ax = PathFinder.Axis(fig[1,1], 
            xlabel = "number of quadrature points per contour (sqrt scale)",
            ylabel = "absolute error (log scale)",
            yscale = log10, xscale = sqrt)

for (i,a) in enumerate(avalsN)
    PathFinder.scatterlines!(nquads, e[i,:], label = "a = $a")
end
PathFinder.axislegend(ax)
PathFinder.limits!(2,26,1e-15,1e-1)
fig