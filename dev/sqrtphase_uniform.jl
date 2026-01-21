using PathFinder
using QuadGK
using CairoMakie
using SpecialFunctions

quadtype = :gaussian
quads = [10,15,30]
# ωvals = [1,5,10,20]
ω = 50

singular = false

avals = 10 .^ range(-5, 0, length=50)
bfix  = 0.0
relErrA = zeros(length(avals), length(quads))

bvals = range(-1.0,1.0, length = 50)
afix  = 1.0
relErrB = zeros(length(bvals), length(quads))

# for ω in ωvals
    
    for (j,a) in enumerate(avals)

        fA(z) = singular ? hankelh1(0, ω*sqrt(z^2 + a^2)) * cis(-ω * sqrt(z^2 + a^2)) : 1.0
        ζA(z) = fA(z) * cis(ω * (sqrt(z^2+a^2) + bfix*z))
        refA = quadgk(ζA, 0.0, 1.0, atol = 1e-14)[1] # brute-force
        G = SquareRootPhaseFunction(a,bfix)
        for (n,N) in enumerate(quads)
            intA = integrate(0.0, 1.0, fA, G, ω; N)[1]   
            relErrA[j,n] = abs(intA-refA)/ abs(refA)      
        end   
    end
    
    
    for (j,b) in enumerate(bvals)
        fB(z) = singular ? hankelh1(0, ω*sqrt(z^2 + afix^2)) * cis(-ω * sqrt(z^2 + afix^2)) : 1.0
        ζB(z) = fB(z) * cis(ω * (sqrt(z^2+afix^2) + b*z))
        refB = quadgk(ζB, 0.0, 1.0, atol = 1e-14)[1] # brute-force
        G = SquareRootPhaseFunction(afix,b)
        for (n,N) in enumerate(quads)
            intB = integrate(0.0, 1.0, fB, G, ω; N)[1]   
            relErrB[j,n] = abs(intB-refB)/ abs(refB)      
        end   
    end
          
# end

fig = Figure(size = (800,400))
axA = Axis(fig[1, 1], title = "Relative error for a ∈ (0,1) with b = $bfix",  
                xlabel = "a (log scale)", ylabel = "Rel. error", limits = (nothing, nothing,1e-16,1e-1), yscale = log10, xscale = log10)
axB = Axis(fig[1, 2], title = "Relative error for b ∈ [-1,1] with a = $afix",  
                xlabel = "b", ylabel = "Rel. error", limits = (-1,1,1e-16,1e-1), yscale = log10)
for (i,n) in enumerate(quads)
    lines!(axA, avals, relErrA[:,i], label = "n = $n")
    lines!(axB, bvals, relErrB[:,i], label = "n = $n")
end
axislegend(axA)
axislegend(axB)
fig