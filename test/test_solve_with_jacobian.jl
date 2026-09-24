# `solve_with_jacobian`, `RCSolution.jacobian`, and restarting `solve` from an `RCSolution`.
using Test, ForwardDiff, StaticArrays
using RheologyCalculator.RheologyModels
import RheologyCalculator: SeriesModel, initial_guess_x, normalisation_x, solve,
    compute_residual, second_invariant_value, _direct_leaf_elastic_correction
import RheologyCalculator.RheologyModels: DruckerPragerCap

@testset "solve: warm start and the converged Jacobian" begin

    c = SeriesModel(
        PowerLawViscosity(1.0e23, 3.0), Elasticity(1.0e10, 2.0e11),
        DruckerPragerCap(; C = 1.0e6, ϕ = 30.0, ψ = 10.0, η_vp = 1.0e19, Pt = -5.0e5),
    )
    others = (; dt = 1.0e8, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,))

    @testset "a solve restarted from its solution converges in fewer iterations" begin
        for rate in (1.0e-13, 1.0e-12)
            vars = vars_2D(rate, 0.0)
            xnorm = normalisation_x(c, 1.0e6, second_invariant_2D(vars.ε) + abs(vars.θ))

            cold = solve(c, initial_guess_x(c, vars, (;), others), vars, others; xnorm0 = xnorm)
            warm = solve(c, cold, vars, others; xnorm0 = xnorm)

            @test warm.iterations < cold.iterations
            @test warm.x ≈ cold.x rtol = 1.0e-12
        end
    end

    @testset "solve_with_jacobian returns the Jacobian at the converged x" begin
        for rate in (1.0e-13, 1.0e-12, 1.0e-11)
            vars = vars_2D(rate, 0.0)
            xnorm = normalisation_x(c, 1.0e6, second_invariant_2D(vars.ε) + abs(vars.θ))
            x0 = initial_guess_x(c, vars, (;), others)
            sol = solve_with_jacobian(c, x0, vars, others; xnorm0 = xnorm)

            # `solve` iterates on the corrected invariant of ε
            ε_corr = _direct_leaf_elastic_correction(c, vars.ε, others)
            v = (; ε = second_invariant_value(vars.ε .+ ε_corr), θ = vars.θ)
            J = ForwardDiff.jacobian(y -> compute_residual(c, y, v, others), sol.x)

            @test sol.jacobian !== nothing
            @test size(sol.jacobian) == size(J)
            @test sol.jacobian ≈ J rtol = 1.0e-12

            # Adding the Jacobian must not disturb the iteration itself.
            plain = solve(c, x0, vars, others; xnorm0 = xnorm)
            @test plain.jacobian === nothing
            @test sol.x == plain.x
            @test sol.iterations == plain.iterations
            @test sol.residual == plain.residual

            # Accepting an RCSolution as the starting point works too.
            @test solve_with_jacobian(c, plain, vars, others; xnorm0 = xnorm).x ≈ sol.x
        end
    end

    @testset "the vector interface is unchanged" begin
        vars = vars_2D(1.0e-12, 0.0)
        xnorm = normalisation_x(c, 1.0e6, second_invariant_2D(vars.ε))
        sol = solve(c, initial_guess_x(c, vars, (;), others), vars, others; xnorm0 = xnorm)
        @test sol isa AbstractVector
        @test collect(sol) == collect(sol.x)
        @test sol[1] == sol.x[1]
        @test sol.iterations isa Int
        @test solve(c, sol, vars, others; xnorm0 = xnorm) ≈ sol
    end
end
