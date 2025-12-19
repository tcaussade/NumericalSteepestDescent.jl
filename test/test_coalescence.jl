using Test
using QuadGK

function coalescence_test(numQuadPts, showText=false)
    maxErr = 0.0
    a = -1.0
    b = 1.0
    freqRange = [10.0, 100.0, 1000.0]
    rArray = 10 .^ range(-4, 0, length=20)
    for freq in freqRange
        if showText
            println("freq=$freq")
        end
        for P in [2,4,6]
            if showText
                println("\tdegree=$P")
            end
            for r in rArray
                if showText
                    println("\t\tr=$r")
                end
                poly_coeffs = [1/(P+1); zeros(P-1); -r^P; 0]
                G = PolynomialPhaseFunction(poly_coeffs)
                f = x -> 1.0
                I_PF = integrate(a, b, f, G, freq; N=numQuadPts)[1]
                integrand(x) = exp(1im * freq * Polynomials.polyval(Polynomials.Polynomial(poly_coeffs), x))
                I_ML = quadgk(integrand, a, b)[1]
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