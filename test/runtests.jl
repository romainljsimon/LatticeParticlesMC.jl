using LatticeParticlesMC
using Test
using StaticArrays

@testset "Distances" begin
    x = SVector{3, Int}(0, 0, 1)
    y1 = SVector{3, Int}(0, 0, 31)
    y2 = SVector{3, Int}(105, 100, 1)
    box = SVector{3, Int}(10, 10, 10)
    r1 = nearest_image_distance_squared(x, y1, box)
    r2 = nearest_image_distance_squared(x, y2, box)
    @test r1 == 0
    @test r2 == 25
    # Write your tests here.
end
