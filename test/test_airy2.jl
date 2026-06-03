using Test
using SpecialFunctions
using NumericalSteepestDescent

"""
    Test the Airy function implementation against Julia's SpecialFunctions.airyai
    See Integral representation in DLMF 9.5.6
"""

function airy_test(nQuadPts, min_z, outputText=false)
    nxPts = 100
    zArray = range(min_z, 5, length=nxPts)
    if outputText
        println("Testing against Julia's Airy function")
    end
    maxRelErr = 0.0
    for z in zArray
        if outputText
            println("z=$z")
        end
        airyJulia = airyai(z)
        G = RationalPhase([0,0,0,im/3],[0],[[0,0,im*z^3/3]])
        integral = nsd([1e-2, 0.0], z -> 1.0, G, 1.0; N=nQuadPts, infcontour=[false,true])
        airyNumerical = (sqrt(3)/(2*π)) * integral
        relErr = abs(airyJulia - airyNumerical) / abs(airyJulia)
        if outputText
            println("\trel err=$relErr")
        end
        maxRelErr = max(maxRelErr, relErr)
    end
    return maxRelErr
end

@testset "Airy Test" begin
    @test airy_test(15, 1) < 1e-6
    @test airy_test(50, 1) < 1e-12
    # For small z, stationary points are very close to pole
    @test airy_test(15, 0.1) < 1e-2
    @test airy_test(50, 0.1) < 1e-3
end

# airy_test(50, 0.1)

# "To debug!"
# z = 0.1
# G = RationalPhase([0,0,0,im/3],[0],[[0,0,im*z^3/3]])
# airyai(z)
# _,fig=integrate(1e-2, 0.0, x -> 1.0, G, 1.0; infcontour=[false,true], plot_sd = true) #; infcontour=[false,true])