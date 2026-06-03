using Test
using SpecialFunctions

"""
    Test the Airy function implementation against Julia's SpecialFunctions.airyai
    See Section 5.2 of the NumericalSteepestDescent paper for details.
"""

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
        G = PolynomialPhase(-1im * [0, -x, 0, 1/3])
        integral = nsd([-π/3, π/3], x -> 1.0, G, 1.0; N=nQuadPts, infcontour=[true,true])
        airyNumerical = (1/(2im*π)) * integral
        relErr = abs(airyJulia - airyNumerical) / abs(airyJulia)
        if outputText
            println("\trel err=$relErr")
        end
        maxRelErr = max(maxRelErr, relErr)
    end
    return maxRelErr
end

@testset "Airy Test" begin
    @test airy_test(15) < 1e-5
    @test airy_test(50) < 1e-10
    @test airy_test(101) < 1e-12
end