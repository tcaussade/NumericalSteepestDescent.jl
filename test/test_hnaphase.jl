using Test
using QuadGK

"""
    Test the SquareRootPhaseFunction struct and evaluation for this phase
"""

function hna_test(nQuadPts, atol, quadtype = :gaussian, outputText=false)
    nvals = 20
    avals = 10 .^ range(-4, 0, length=nvals)
    bvals = range(-1.0,1.0, length = nvals)
    ωvals = [1,5,10,20]
    if outputText
        println("Testing self-error for SquareRoot Phase")
    end
    maxRelErr = 0.0

    for ω in ωvals
        for a in avals
            for b in bvals
                ζ(z) = cis(ω * (sqrt(z^2+a^2) + b*z))
                if outputText
                    println("ω = $ω, a = $a, b = $b")
                end
                ref = quadgk(ζ, 0.0, 1.0, atol = 1e-14)[1] # brute-force
                G = SquareRootPhaseFunction(a,b)
                if quadtype == :gaussian
                    int = integrate(0.0, 1.0, x -> 1.0, G, ω; N=nQuadPts)[1]            
                elseif quadtype == :adaptive
                    int = integrate(0.0, 1.0, x -> 1.0, G, ω; quadtype, atol)[1]
                end

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

@testset "HNA Phase Test" begin
    @testset "Gaussian quadrature" begin
        @test hna_test(15,0, :gaussian) < 1e-3
        @test hna_test(50,0, :gaussian) < 1e-5
        @test hna_test(101,0, :gaussian) < 3e-7
    end
    # @testset "Adaptive quadrature" begin
    #     @test hna_test(0, 1e-2, :adaptive) < 1e-2
    #     @test hna_test(0, 1e-6, :adaptive) < 1e-5
    #     @test hna_test(0, 1e-10, :adaptive) < 1e-10
    # end
end

# hna_test(0, 1e-12, :adaptive) 
# hna_test(100, 0, :gaussian, false)