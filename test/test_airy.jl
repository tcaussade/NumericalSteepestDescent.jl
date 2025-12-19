using Test
using SpecialFunctions

function airy_test(nQuadPts, outputText=false)
    nxPts = 100
    xArray = range(-5, 5, length=nxPts)
    if outputText
        println("Testing against Julia's Airy function")
    end
    maxRelErr = 0.0
    for x in xArray
        if outputText
            println("x=$x")
        end
        airyJulia = airyai(x)
        coeffs = -1im * [1/3, 0, -x, 0]
        G = PolynomialPhaseFunction(coeffs)
        f = x -> 1.0
        ω = 1.0
        a = -π/3
        b = π/3
        integral = integrate(a, b, f, G, ω; N=nQuadPts)[1]
        airyPathFinder = (1/(2im*π)) * integral
        relErr = abs(airyJulia - airyPathFinder) / abs(airyJulia)
        if outputText
            println("\trel err=$relErr")
        end
        maxRelErr = max(maxRelErr, relErr)
    end
    return maxRelErr
end

@testset "Airy Test" begin
    @test airy_test(15) < 1e-6
    @test airy_test(50) < 1e-10
    @test airy_test(101) < 1e-10
end