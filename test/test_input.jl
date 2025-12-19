using Test

function test_input_check()
    isPass = true
    # Bad a
    try
        integrate([0, 0], 1, x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch
    end
    # Bad b
    try
        integrate(-1, [0, 0], x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch
    end
    # Bad f
    try
        integrate(-1, 1, "oops", PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10)
        isPass = false
    catch
    end
    # Bad coeffs
    try
        integrate(-1, 1, x -> x^2, PolynomialPhaseFunction(rand(2)), 50; N=10)
        isPass = false
    catch
    end
    # Bad ω
    try
        integrate(-1, 1, x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), [1, 1]; N=10)
        isPass = false
    catch
    end
    # Bad N
    try
        integrate(-1, 1, x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=-10)
        isPass = false
    catch
    end
    return isPass
end

@testset "Input Check" begin
    @test test_input_check()
end