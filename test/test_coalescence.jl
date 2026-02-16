using Test
using QuadGK

"""
    See Section 5.4 of the PathFinder paper for details.
"""

function coalescence_test(numQuadPts, showText=false)
    maxErr = 0.0
    a,b = (-1.0, 1.0)
    freqRange = [10.0, 100.0, 1000.0]
    rArray = 10 .^ range(-4, 0, length=20) # when r = 0 there is coalescence
    for freq in freqRange
        if showText
            println("freq=$freq")
        end
        for P in [2,4,6] # order of stationary point
            if showText
                println("\tdegree=$P")
            end
            for r in rArray
                if showText
                    println("\t\tr=$r")
                end
                poly_coeffs = [0.0; -r^P; zeros(P-1); 1.0/(P+1)]  
                G = PolynomialPhaseFunction(poly_coeffs)
                g(z) = PathFinder.evalphase(z,G)
                I_PF,_ = integrate(a, b, z -> 1.0, G, freq; N=numQuadPts) 
                I_ML,_ = quadgk(x -> cis(freq * g(x)), a, b)
                Ierr = abs(I_PF - I_ML) / abs(I_ML)
                if showText
                    println("\t\trel err=$Ierr")
                end
                maxErr = max(maxErr, Ierr)
            end
        end
    end
    return maxErr
end

@testset "Coalescence Test" begin
    @test coalescence_test(25) < 1e-4
    @test coalescence_test(75) < 1e-12
    @test coalescence_test(101) < 1e-12
end