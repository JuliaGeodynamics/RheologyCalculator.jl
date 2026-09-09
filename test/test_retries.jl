using Test
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "host-side retries" begin
    c = SeriesModel(LinearViscosity(1.0e22))
    vars = (; ε = 1.0e-14)
    others = (;)
    x = initial_guess_x(c, vars, (; τ = 0.0), others)

    @test solve_with_retries(c, x, vars, others; max_retries = 0).x ==
        solve(c, x, vars, others).x
    @test_throws ArgumentError solve_with_retries(c, x, vars, others; max_retries = -1)
    @test_throws ArgumentError solve_with_retries(c, x, vars, others; tolerance_factor = 0.5)
end
