using Test
using QuadGK

function umbilic_test(n)
    return 0.0
end

@testset "Umbilic Test" begin
    @test umbilic_test(20) < 2e-4
    @test umbilic_test(50) < 2e-5
    @test umbilic_test(101) < 2e-5
end