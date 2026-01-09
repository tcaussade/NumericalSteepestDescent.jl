using Test
using SpecialFunctions

"""
    Test the SquareRootPhaseFunction struct and evaluation for this phase
"""

function hna_test(nQuadPts, outputText=false)
    # nxPts = 100
    # xArray = range(-5, 5, length=nxPts)
    # if outputText
    #     println("Testing against Julia's Airy function")
    # end
    # maxRelErr = 0.0
    # for x in xArray
    #     if outputText
    #         println("x=$x")
    #     end
    #     airyJulia = airyai(x)
    #     G = PolynomialPhaseFunction(-1im * [0, -x, 0, 1/3])
    #     integral,_ = integrate(-π/3, π/3, x -> 1.0, G, 1.0; N=nQuadPts, infcontour=[true,true])
    #     airyPathFinder = (1/(2im*π)) * integral
    #     relErr = abs(airyJulia - airyPathFinder) / abs(airyJulia)
    #     if outputText
    #         println("\trel err=$relErr")
    #     end
    #     maxRelErr = max(maxRelErr, relErr)
    # end
    # return maxRelErr
end

@testset "HNA Phase Test" begin
    @test hna_test(15) < 1e-6
    @test hna_test(50) < 1e-10
    @test hna_test(101) < 1e-10
end