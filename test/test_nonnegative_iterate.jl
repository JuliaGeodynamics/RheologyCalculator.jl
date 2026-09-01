using RheologyCalculator, Test, StaticArrays
using RheologyCalculator.RheologyModels
import RheologyCalculator.RheologyModels: second_invariant_2D
import RheologyCalculator: max_feasible_step, branch_strain_rate_mask

@testset "non-negative branch strain-rate iterate" begin
    # A power law in parallel is undefined for a negative branch strain rate,
    # and the unbounded Newton iterate overshoots below zero at low imposed
    # strain rate. The step bound keeps the iterate inside the physical range
    # across the whole range, not just where the overshoot happens to be small.
    c = SeriesModel(
        LinearViscosity(1.0e3),
        ParallelModel(PowerLawViscosity(4.0e-4, 3), LinearViscosity(1.0)),
    )

    function run(e)
        εᵢⱼ = (e, -e, 0.0)
        εII = second_invariant_2D(εᵢⱼ)
        g = initial_guess_x(c, (; ε = εII, θ = 0.0), (; τ = 0.0), (;))
        return solve(c, g, (; ε = εᵢⱼ, θ = 0.0), (;); xnorm0 = normalisation_x(c, 1.0e6, εII))
    end

    for e in (1.0e-3, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-10, 1.0e-15, 1.0e-20)
        x = run(e)
        @test all(isfinite, x)
        @test all(≥(0), x)
    end

    # Well below the elastic/power-law crossover the linear element carries the
    # whole strain rate, so τII → 2ηεII and the branch strain rate vanishes.
    x = run(1.0e-20)
    εII = second_invariant_2D((1.0e-20, -1.0e-20, 0.0))
    @test x[1] ≈ 2 * 1.0e3 * εII rtol = 1.0e-6

    # The strain rates of the two sides of the branch sum to the imposed one.
    x = run(1.0e-3)
    εII = second_invariant_2D((1.0e-3, -1.0e-3, 0.0))
    @test x[1] / (2 * 1.0e3) + x[2] ≈ εII rtol = 1.0e-10
end

@testset "max_feasible_step" begin
    mask = SA[true, false]

    # No live constraint: the full step is allowed.
    @test max_feasible_step(SA[1.0, 1.0], SA[1.0, -10.0], mask) == 1.0
    @test max_feasible_step(SA[1.0, 1.0], SA[-0.1, -10.0], mask) == 1.0

    # A step that would cross zero is cut back to just short of the boundary.
    α = max_feasible_step(SA[1.0, 1.0], SA[-2.0, 0.0], mask)
    @test α ≈ 0.995 / 2
    @test (SA[1.0, 1.0] + α * SA[-2.0, 0.0])[1] > 0

    # An unmasked entry is never bounded, however far negative it goes.
    @test max_feasible_step(SA[1.0, 1.0], SA[0.0, -100.0], mask) == 1.0

    # An entry already at the boundary imposes no bound, so the step is never
    # cut to zero and the iteration cannot stall here.
    @test max_feasible_step(SA[0.0, 1.0], SA[-1.0, 0.0], mask) == 1.0
end

@testset "branch_strain_rate_mask" begin
    # Only a parallel branch's strain-rate unknown is marked.
    c = SeriesModel(
        LinearViscosity(1.0e3),
        ParallelModel(PowerLawViscosity(4.0e-4, 3), LinearViscosity(1.0)),
    )
    @test x_keys(c) === (:τ, :ε)
    @test branch_strain_rate_mask(c) == SA[false, true]

    # A composite with no parallel branch has no bounded entry.
    c_series = SeriesModel(LinearViscosity(1.0e22), Elasticity(1.0e10, 1.0e11))
    @test !any(branch_strain_rate_mask(c_series))
    @test length(branch_strain_rate_mask(c_series)) == length(generate_equations(c_series))
end
