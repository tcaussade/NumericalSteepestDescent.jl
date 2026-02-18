using Test
using QuadGK

"""
    Test the SquareRootPhaseFunction struct and evaluation for this phase
"""


function testparameters(outputText)
    nvals = 10
    avals = 10 .^ range(-8, 0, length=nvals)
    bvals = range(-1.0,1.0, length = nvals)
    ωvals = [1,5,10,20]
    if outputText
        println("Testing self-error for SquareRoot Phase")
    end
    return ωvals, avals, bvals
end

function hna_test_gaussian(nQuadPts, outputText=false)
    ωvals, avals, bvals = testparameters(outputText)
    maxAbsErr = 0.0
    for ω in ωvals
        for a in avals
            for b in bvals[2:end-1] # remove b=±1
                ζ(z) = cis(ω * (sqrt(z^2+a^2) + b*z))
                ref = quadgk(ζ, 0.0, 1.0, atol = 1e-14)[1] # brute-force
                if outputText
                    println("ω = $ω, a = $a, b = $b")
                end          
                G = SquareRootPhaseFunction(a,b)
                int = integrate(0.0, 1.0, x -> 1.0, G, ω; N=nQuadPts)            
                absErr = abs(int-ref) #/ abs(ref)
                if outputText
                    println("\t abs err = $absErr")
                end
                maxAbsErr = max(maxAbsErr, absErr)
            end
        end
    end
    return maxAbsErr
end

function hna_test_adaptive(atol, outputText=false)
    avals, bvals = testparameters()
    maxRelErr = 0.0
    for ω in ωvals
        for a in avals
            for b in bvals[2:end-1] # remove b=±1
                ζ(z) = cis(ω * (sqrt(z^2+a^2) + b*z))
                ref = quadgk(ζ, 0.0, 1.0, atol = 1e-14)[1] # brute-force
                if outputText
                    println("ω = $ω, a = $a, b = $b")
                end
                G = SquareRootPhaseFunction(a,b)          
                int = integrate(0.0, 1.0, x -> 1.0, G, ω; quadtype, atol)[1]
                relErr = abs(int-ref)/ abs(ref)
                if outputText
                    println("\t rel err = $relErr")
                end
                maxRelErr = max(maxRelErr, relErr)
            end
        end
    end
    return maxRelErr
end


@testset "Gaussian quadrature" begin
    @test hna_test_gaussian(25) < 1e-2
    @test hna_test_gaussian(50) < 1e-5
    @test hna_test_gaussian(100) < 3e-10
end
@testset "Adaptive quadrature" begin
#     @test hna_test_adaptive(1e-2) < 1e-2
#     @test hna_test_adaptive(1e-6) < 1e-5
#     @test hna_test_adaptive(1e-10) < 1e-10
end


# hna_test(0, 1e-12, :adaptive) 