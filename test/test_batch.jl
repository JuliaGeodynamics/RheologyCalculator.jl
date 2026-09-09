using Test
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "static batch solve" begin
    c = SeriesModel(LinearViscosity(1.0e22))
    vars = ((; ε = 1.0e-14), (; ε = 2.0e-14))
    others = ((;), (;))
    xs = (
        initial_guess_x(c, vars[1], (; τ = 0.0), others[1]),
        initial_guess_x(c, vars[2], (; τ = 0.0), others[2]),
    )

    batch = solve_batch(c, xs, vars, others)
    scalar = (solve(c, xs[1], vars[1], others[1]), solve(c, xs[2], vars[2], others[2]))

    @test batch == scalar
    @test isbitstype(typeof(batch))
    solve_batch(c, xs, vars, others)
    @test @allocated(solve_batch(c, xs, vars, others)) == 0
end
