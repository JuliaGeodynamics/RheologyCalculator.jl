# `bt_line_search`'s `(α, x_next, r)` must stay internally consistent.
using Test, StaticArrays
using RheologyCalculator.RheologyModels
import RheologyCalculator: SeriesModel, ParallelModel, initial_guess_x,
    normalisation_x, solve, compute_residual, mynorm, bt_line_search, jacobian,
    backsolve, second_invariant_value, _direct_leaf_elastic_correction

@testset "line search residual reuse" begin

    c = SeriesModel(
        LinearViscosity(1.0e22), IncompressibleElasticity(10.0e9),
        ParallelModel(LinearViscosity(1.0e20), DruckerPrager(10.0e6, 30.0, 0.0)),
    )
    vars = vars_2D(1.0e-14)
    others = (; dt = 1.0e8, τ0 = (zero_stress_tensor_2D(),), P = 1.0e6, P0 = (0.0,))
    xnorm = normalisation_x(c, 1.0e6, second_invariant_2D(vars.ε))
    x0 = initial_guess_x(c, vars, (; τ = 2.0e3, λ = 0.0), others)

    # `solve` iterates on the corrected invariant of ε
    ε_corr = _direct_leaf_elastic_correction(c, vars.ε, others)
    v = (; ε = second_invariant_value(vars.ε .+ ε_corr), θ = vars.θ)

    @testset "the returned triple is self-consistent" begin
        x = x0
        for _ in 1:5
            r = compute_residual(c, x, v, others)
            er = mynorm(r, xnorm)
            Δx = backsolve(jacobian(c, x, v, others), r)

            α, x_next, r_next = bt_line_search(
                Δx, x, c, v, others, xnorm, er; α = 1.0, ρ = 0.5, lstol = 0.95, α_min = 0.1
            )

            @test x_next ≈ x .+ α .* Δx rtol = 1.0e-14
            @test r_next ≈ compute_residual(c, x_next, v, others) rtol = 1.0e-14
            # while the update is significant the step must move; once
            # converged, Δx ≈ 0 and x_next == x legitimately
            if maximum(abs, Δx) > 8 * eps() * maximum(abs, x)
                @test x_next != x
            end
            x = x_next
        end
    end

    @testset "the solve still converges on the hard step" begin
        x = x0
        τ_e, P_e = (zero_stress_tensor_2D(),), (0.0,)
        x_hard, o_hard = x, others
        for _ in 1:469
            o = (; dt = others.dt, τ0 = τ_e, P = others.P, P0 = P_e)
            x_hard, o_hard = x, o
            x = solve(c, x, vars, o; xnorm0 = xnorm)
            τ_e = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, o)
        end
        sol = solve(c, x_hard, vars, o_hard; xnorm0 = xnorm, itermax = 10)
        @test sol.residual < 1.0e-8
    end
end
