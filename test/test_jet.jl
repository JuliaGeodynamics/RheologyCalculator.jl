using JET
import RheologyCalculator: compute_pressure_elastic, compute_residual, compute_stress_elastic
include("type_stability_fixtures.jl")

# Run with `test_args = ["--jet-only"]`: JET only supports released Julia
# versions, so the normal suite leaves it out.
@testset "JET" begin
    for (name, c, vars, args, others) in type_stability_fixtures()
        @testset "$name" begin
            x = initial_guess_x(c, vars, args, others)
            # `target_modules` scopes the report to our own code, ignoring
            # dynamic dispatch inside Base's type-printing / StaticArrays paths.
            JET.@test_opt target_modules = (RheologyCalculator,) generate_equations(c)
            JET.@test_opt target_modules = (RheologyCalculator,) initial_guess_x(c, vars, args, others)
            JET.@test_opt target_modules = (RheologyCalculator,) normalisation_x(c, 1.0e6, 1.0e-15)
            JET.@test_opt target_modules = (RheologyCalculator,) compute_residual(c, x, vars, others)
            JET.@test_opt target_modules = (RheologyCalculator,) compute_stress_elastic(c, x, others)
            JET.@test_opt target_modules = (RheologyCalculator,) compute_pressure_elastic(c, x, others)
        end
    end
end
