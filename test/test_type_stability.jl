using JET
import RheologyCalculator: compute_pressure_elastic, compute_residual, compute_stress_elastic

function type_stability_fixture()
    viscous = LinearViscosity(5.0e19)
    elastic = Elasticity(1.0e10, 1.0e12)
    c = SeriesModel(elastic, viscous)

    vars = (; ε = 1.0e-15, θ = 1.0e-20)
    args = (; τ = 1.0e2, P = 1.0e6)
    others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    x = initial_guess_x(c, vars, args, others)
    xnorm = normalisation_x(c, 1.0e6, 1.0e-15)

    return c, vars, others, x, xnorm
end

@testset "Type stability" begin
    c, vars, others, x, xnorm = type_stability_fixture()

    eqs = @inferred generate_equations(c)
    @test length(eqs) == length(x)

    @inferred initial_guess_x(c, vars, (; τ = 1.0e2, P = 1.0e6), others)
    @inferred normalisation_x(c, 1.0e6, 1.0e-15)
    @inferred compute_residual(c, x, vars, others)
    @inferred compute_stress_elastic(c, x, others)
    @inferred compute_pressure_elastic(c, x, others)
    @inferred solve(c, x, vars, others; xnorm0 = xnorm, itermax = 1, verbose = false)

    # JET's optimization analysis is sensitive to compiler-internal inference
    # changes and is unreliable on pre-release Julia, so only run it on released
    # versions. `target_modules` scopes the report to our own code, ignoring
    # dynamic dispatch inside Base's type-printing / StaticArrays paths.
    if isempty(VERSION.prerelease)
        JET.@test_opt target_modules = (RheologyCalculator,) generate_equations(c)
        JET.@test_opt target_modules = (RheologyCalculator,) initial_guess_x(c, vars, (; τ = 1.0e2, P = 1.0e6), others)
        JET.@test_opt target_modules = (RheologyCalculator,) normalisation_x(c, 1.0e6, 1.0e-15)
        JET.@test_opt target_modules = (RheologyCalculator,) compute_residual(c, x, vars, others)
        JET.@test_opt target_modules = (RheologyCalculator,) compute_stress_elastic(c, x, others)
        JET.@test_opt target_modules = (RheologyCalculator,) compute_pressure_elastic(c, x, others)
    end
end
