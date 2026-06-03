using Test
using NumericalSteepestDescent
using SpecialFunctions
using QuadGK

@testset "NumericalSteepestDescent Tests" begin
    include("test_plot.jl")
    include("test_input.jl")
    @testset "Polynomial phase" begin
        include("test_airy.jl")
        include("test_coalescence.jl")
        include("test_pearcey.jl")
    end
    @testset "Square-Root phase" begin
        include("test_hnaphase.jl")
    end
    @testset "Rational phase" begin
        include("test_airy2.jl")
        # include("test_umbilic.jl")
    end
end