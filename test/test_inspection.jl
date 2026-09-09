using Test
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "inspection metadata" begin
    c = SeriesModel(
        LinearViscosity(1.0e22),
        ParallelModel(LinearViscosity(1.0e21), IncompressibleElasticity(1.0e10)),
    )
    entries = inspect(c)

    @test entries[1].parent == 0
    @test entries[1].children == (2,)
    @test entries[1].inputs == (:ε,)
    @test entries[1].history == (:d,)
    @test entries[2].inputs == (:τ,)
    @test entries[2].history == (:d, :τ0, :P0)
end
