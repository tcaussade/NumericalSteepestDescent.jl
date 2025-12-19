using Test
using PathFinder
using SpecialFunctions
using QuadGK

@testset "PathFinder Tests" begin
    include("test_airy.jl")
    include("test_coalescence.jl")
    include("test_pearcey.jl")
    include("test_plot.jl")
    include("test_input.jl")
end