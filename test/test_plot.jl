using Test

function plot_test()
    try
        integrate(-1, 1, x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10, plot_graph=true)
        integrate(-1, 1, x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10, plot_sd=true)
        integrate(-1, 1, x -> x^2, PolynomialPhaseFunction([1, -0.5, 0.5, 0, -1, 0]), 50; N=10, plot_graph=true, plot_sd=true)
        integrate(-1, 1, x -> x^2, PolynomialPhaseFunction([1, 0]), 50; N=10, plot_graph=true, plot_sd=true)
        return true
    catch
        return false
    end
end

@testset "Plot Test" begin
    @test plot_test()
end