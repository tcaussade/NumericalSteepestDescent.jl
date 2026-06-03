using Test

function test_input_check_polynomials()
    isPass = true  

    try # Bad endpoints
        nsd([0, 0], 1, x -> x^2, PolynomialPhase([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    
    try 
        nsd(-1, [0, 0], x -> x^2, PolynomialPhase([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    
    try # Bad f
        nsd(-1, 1, "oops", PolynomialPhase([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end

    try 
        nsd(-1, 0, x -> x^2, PolynomialPhase([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    
    try # Bad coeffs
        nsd([-1, 1], x -> x^2, PolynomialPhase(rand(2)), 50; N=10)
        isPass = false
    catch end
    
    try # Bad ω
        nsd([-1, 1], x -> x^2, PolynomialPhase([1, -0.5, 0.5, 0, -1, 0]), [1, 1]; N=10)
        isPass = false
    catch end
    
    try # Bad N
        nsd([-1, 1], x -> x^2, PolynomialPhase([1, -0.5, 0.5, 0, -1, 0]), 50; N=-10)
        isPass = false
    catch end

    return isPass
end

function test_input_check_linear()
    isPass = true
    try # Bad calling
        nsd([-1, 1], x -> x^2, LinearPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    return isPass
end

function test_input_check_sqrt()
    isPass = true
    try # Bad calling
        nsd([-1, 1], x -> x^2, SquareRootPhase([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    return isPass
end


@testset "Input Check" begin
    @testset "PolynomialPhase" begin
        @test test_input_check_polynomials()
    end
    @testset "LinearPhaseFunction" begin 
        @test test_input_check_linear()
    end
    @testset "SquareRootPhase" begin
        @test test_input_check_sqrt()
    end
end