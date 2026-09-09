using Test
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "GPU-facing invariants" begin
    c = SeriesModel(
        LinearViscosity(1.0e22),
        ParallelModel(LinearViscosity(1.0e21), IncompressibleElasticity(1.0e10)),
    )
    vars = (; ε = 1.0e-14, θ = 0.0)
    others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    x = initial_guess_x(c, vars, (; τ = 0.0, P = 0.0), others)
    sol = solve(c, x, vars, others)

    @test isbitstype(typeof(c))
    @test isbitstype(typeof(vars))
    @test isbitstype(typeof(others))
    @test isbitstype(typeof(x))
    @test isbitstype(typeof(sol))
    @test isbitstype(typeof(jacobian(c, x, vars, others)))
end
