using RheologyCalculator, Test
using RheologyCalculator.RheologyModels

@testset "solve with throw_on_failure = false" begin
    c = SeriesModel(LinearViscosity(1.0e22), ParallelModel(LinearViscosity(1.0e21), IncompressibleElasticity(1.0e10)))
    vars = (; ε = 1.0e-14)
    others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    x0 = initial_guess_x(c, vars, (; τ = 0.0), others)

    @testset "converged solve is unchanged" begin
        sol = solve(c, x0, vars, others)
        sol_nothrow = solve(c, x0, vars, others; throw_on_failure = false)
        @test sol.converged
        @test sol_nothrow.converged
        @test sol_nothrow.x == sol.x
        @test sol_nothrow.iterations == sol.iterations
        @test solve_with_jacobian(c, x0, vars, others).converged
    end

    @testset "iteration limit" begin
        rs = RheologyModels.RateStateFriction(0.0, 0.2, 4.0e-9, 0.011, 0.017, 0.05, 0.0, 625.0)
        cr = SeriesModel(IncompressibleElasticity(3.0e10), rs)
        vr = (; ε = 1.0e-9)
        or = (; dt = 1.0e6, τ0 = (0.0,), P0 = (0.0,), P = 5.0e7, Ω_old = -5.0)
        xr = initial_guess_x(cr, vr, (; τ = 0.0, P = 5.0e7), or)
        kw = (; xnorm0 = normalisation_x(cr, 1.0e6, 1.0e-9), itermax = 2)
        err = try
            solve(cr, xr, vr, or; kw...)
        catch e
            e
        end
        @test err isa NonConvergenceError
        sol = solve(cr, xr, vr, or; kw..., throw_on_failure = false)
        @test !sol.converged
        @test sol.x == err.x
        @test sol.iterations == err.iterations
        @test sol.residual == err.residual
        @test !solve_with_jacobian(cr, xr, vr, or; kw..., throw_on_failure = false).converged
        @test occursin("not converged", sprint(show, MIME("text/plain"), sol))
    end

    @testset "stagnation" begin
        # Warm-started from its own solution without `xnorm0`, this model sits at
        # a residual floor just above `atol`, so the iteration stagnates.
        v1, v2 = LinearViscosity(5.0e19), LinearViscosity(1.0e20)
        cs = SeriesModel(v1, ParallelModel(SeriesModel(v1, v2), v2), ParallelModel(v1, v2))
        vs = (; ε = 1.0e-15)
        warm = solve(cs, initial_guess_x(cs, vs, (; τ = 0.0), (;)), vs, (;))
        @test_throws NonConvergenceError solve(cs, warm, vs, (;))
        sol = solve(cs, warm, vs, (;); throw_on_failure = false)
        @test !sol.converged
        @test sol.iterations < 10
    end
end
