using Test
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "validation" begin
    c = SeriesModel(LinearViscosity(1.0e22), IncompressibleElasticity(1.0e10))
    vars = (; ε = 1.0e-14, θ = 0.0)
    others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))

    @test validate(c, vars, others) === c
    @test isbitstype(typeof(c))

    @test_throws ArgumentError validate(SeriesModel(LinearViscosity(0.0)), vars, (;))
    @test_throws ArgumentError validate(c, (; ε = NaN, θ = 0.0), others)
    @test_throws ArgumentError validate(c, vars, (; dt = 1.0e10, τ0 = (0.0,)))
end
