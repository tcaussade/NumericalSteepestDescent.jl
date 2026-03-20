using Test

function test_input_check_polynomials()
    isPass = true  

    try # Bad endpoints
        integrate([0, 0], 1, x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    
    try 
        integrate(-1, [0, 0], x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    
    try # Bad f
        integrate(-1, 1, "oops", PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end

    try 
        integrate(-1, 0, x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    
    try # Bad coeffs
        integrate([-1, 1], x -> x^2, PolynomialPhaseFunction(rand(2)), 50; N=10)
        isPass = false
    catch end
    
    try # Bad ω
        integrate([-1, 1], x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), [1, 1]; N=10)
        isPass = false
    catch end
    
    try # Bad N
        integrate([-1, 1], x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=-10)
        isPass = false
    catch end

    return isPass
end

function test_input_check_linear()
    isPass = true
    try # Bad calling
        integrate([-1, 1], x -> x^2, LinearPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    return isPass
end

function test_input_check_sqrt()
    isPass = true
    try # Bad calling
        integrate([-1, 1], x -> x^2, SquareRootPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch end
    return isPass
end


@testset "Input Check" begin
    @testset "PolynomialPhaseFunction" begin
        @test test_input_check_polynomials()
    end
    @testset "LinearPhaseFunction" begin 
        @test test_input_check_linear()
    end
    @testset "SquareRootPhaseFunction" begin
        @test test_input_check_sqrt()
    end
end