using JET
import RheologyCalculator: compute_pressure_elastic, compute_residual, compute_stress_elastic

# Each fixture is (name, composite, vars, args, others). The topologies are the
# ones whose equation generation is unrolled per composite type: losing
# inference on any of them turns a static residual evaluation into a dynamic
# one without changing a single number.
function type_stability_fixtures()
    viscous = LinearViscosity(5.0e19)
    viscous2 = LinearViscosity(1.0e20)
    elastic = Elasticity(1.0e10, 1.0e12)
    elasticinc = IncompressibleElasticity(1.0e10)
    bulkvisc = BulkViscosity(1.0e18)

    return (
        (
            "series elastic-viscous",
            SeriesModel(elastic, viscous),
            (; ε = 1.0e-15, θ = 1.0e-20),
            (; τ = 1.0e2, P = 1.0e6),
            (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,)),
        ),
        (
            "Kelvin-Voigt branch",
            SeriesModel(elasticinc, ParallelModel(elasticinc, viscous)),
            (; ε = 1.0e-15),
            (; τ = 1.0e2),
            (; dt = 1.0e10, τ0 = (0.0, 0.0), P0 = (0.0, 0.0)),
        ),
        (
            "generalized Maxwell branch",
            SeriesModel(viscous, ParallelModel(SeriesModel(viscous2, elasticinc), viscous2)),
            (; ε = 1.0e-15),
            (; τ = 1.0e2),
            (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,)),
        ),
        (
            "plastic series",
            SeriesModel(LinearViscosity(1.0e22), elasticinc, DruckerPrager(1.0e6, 30.0, 0.0)),
            (; ε = 1.0e-14),
            (; τ = 1.0e6),
            (; dt = 1.0e10, τ0 = (0.0,), P = 1.0e6, P0 = (0.0,)),
        ),
        (
            "volumetric series",
            SeriesModel(elastic, bulkvisc),
            (; ε = 1.0e-15, θ = 1.0e-20),
            (; τ = 1.0e2, P = 1.0e6),
            (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,)),
        ),
    )
end

@testset "Type stability" begin
    for (name, c, vars, args, others) in type_stability_fixtures()
        @testset "$name" begin
            x = initial_guess_x(c, vars, args, others)
            xnorm = normalisation_x(c, 1.0e6, 1.0e-15)

            eqs = @inferred generate_equations(c)
            @test length(eqs) == length(x)

            @inferred initial_guess_x(c, vars, args, others)
            @inferred normalisation_x(c, 1.0e6, 1.0e-15)
            @inferred compute_residual(c, x, vars, others)
            @inferred compute_stress_elastic(c, x, others)
            @inferred compute_pressure_elastic(c, x, others)
            @inferred solve(c, x, vars, others; xnorm0 = xnorm, itermax = 1, verbose = false)

            # JET's optimization analysis is sensitive to compiler-internal
            # inference changes and is unreliable on pre-release Julia, so only
            # run it on released versions. `target_modules` scopes the report to
            # our own code, ignoring dynamic dispatch inside Base's type-printing
            # / StaticArrays paths.
            if isempty(VERSION.prerelease)
                JET.@test_opt target_modules = (RheologyCalculator,) generate_equations(c)
                JET.@test_opt target_modules = (RheologyCalculator,) initial_guess_x(c, vars, args, others)
                JET.@test_opt target_modules = (RheologyCalculator,) normalisation_x(c, 1.0e6, 1.0e-15)
                JET.@test_opt target_modules = (RheologyCalculator,) compute_residual(c, x, vars, others)
                JET.@test_opt target_modules = (RheologyCalculator,) compute_stress_elastic(c, x, others)
                JET.@test_opt target_modules = (RheologyCalculator,) compute_pressure_elastic(c, x, others)
            end
        end
    end
end
