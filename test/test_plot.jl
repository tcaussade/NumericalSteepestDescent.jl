using Test

function plot_test()
    try
        nsd([-1, 1], x -> x^2, PolynomialPhase([0, -1, 0, 0.5, -0.5, 1]), 50; N=10, plot_graph=true)
        nsd([-1, 1], x -> x^2, PolynomialPhase([0, -1, 0, 0.5, -0.5, 1]), 50; N=10, plot_sd=true)
        nsd([-1, 1], x -> x^2, PolynomialPhase([0, -1, 0, 0.5, -0.5, 1]), 50; N=10, plot_graph=true, plot_sd=true)
        nsd([-1, 1], x -> x^2, LinearPhase(), 50; N=10, plot_sd=true, plot_graph = true)
        nsd([-1, 1], x -> x^2, SquareRootPhase(1.0, -0.6), 50; N=10, plot_sd=true, plot_graph = true)
        nsd([-1, 1], x -> x^2, SquareRootPhase(1.0, +0.6), 50; N=10, plot_sd=true, plot_graph = true)
        return true
    catch
        return false
    end
end

@testset "Plot Test" begin
    @test plot_test()
end