import RheologyCalculator: compute_pressure_elastic, compute_residual, compute_stress_elastic
include("type_stability_fixtures.jl")

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
        end
    end
end
