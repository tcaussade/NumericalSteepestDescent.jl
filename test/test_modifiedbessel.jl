using Test
using SpecialFunctions

"""
    Test the modified Bessel function implementation against Julia's SpecialFunctions.besselk
    See BesselK definition in DLMF 10.32.10

@warn We need singular quadrature as the amplitude is singular at t=0
"""

function bessel_test(nQuadPts, outputText=false)

    ν = 0

    nxPts = 40
    xArray = range(2, 5, length=nxPts)
    if outputText
        println("Testing against Julia's Modified Bessel function")
    end
    maxRelErr = 0.0
    for x in xArray
        if outputText
            println("x=$x")
        end
        besselJulia = besselk(ν, x)
        G = RationalPhaseFunction(im*[0,1],[0],[im*[x^2/4]])
        f(t) = 1/t^(ν+1)
        integra = integrate(1e-6, 0, f, G, 1.0; N=nQuadPts, infcontour=[false,true])
        besselPathFinder = 0.5*(0.5*x)^ν * integral
        relErr = abs(besselJulia - besselPathFinder) / abs(besselJulia)
        if outputText
            println("\trel err=$relErr")
        end
        maxRelErr = max(maxRelErr, relErr)
    end
    return maxRelErr
end

@testset "Bessel Test" begin
    @test bessel_test(15) < 1e-2
    # @test bessel_test(50) < 1e-10
    # @test bessel_test(101) < 1e-10
end

# bessel_test(15, true)