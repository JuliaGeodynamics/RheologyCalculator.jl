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

@testset "retries relax the tolerances" begin
    # Warm-started from its own solution without `xnorm0`, this model sits at a
    # residual floor of 7.3e-12: it fails atol = 1e-12 and passes 1e-11.
    v1, v2 = LinearViscosity(5.0e19), LinearViscosity(1.0e20)
    c = SeriesModel(v1, ParallelModel(SeriesModel(v1, v2), v2), ParallelModel(v1, v2))
    vars = (; ε = 1.0e-15)
    warm = solve(c, initial_guess_x(c, vars, (; τ = 0.0), (;)), vars, (;)).x

    @test_throws NonConvergenceError solve(c, warm, vars, (;))
    sol = solve_with_retries(c, warm, vars, (;); max_retries = 1)
    @test sol.converged
    @test sol.x ≈ warm rtol = 1.0e-14
    # The caller's throw_on_failure doesn't stop the retries.
    @test solve_with_retries(c, warm, vars, (;); max_retries = 1, throw_on_failure = false).converged
    err = try
        solve_with_retries(c, warm, vars, (;); max_retries = 1, tolerance_factor = 1)
    catch e
        e
    end
    @test err isa NonConvergenceError
    @test err.reason == :stagnation
end
