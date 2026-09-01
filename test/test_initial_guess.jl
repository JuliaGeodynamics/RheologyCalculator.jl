using RheologyCalculator, Test, StaticArrays
using RheologyCalculator.RheologyModels
import RheologyCalculator.RheologyModels: second_invariant_2D

@testset "initial guess for parallel branch strain rates" begin
    # A ParallelModel sibling of another leaf introduces a branch strain-rate
    # unknown. Seeding it at zero is a singular point of the parallel effective
    # viscosity of a power law, which makes the first Newton step NaN.
    c = SeriesModel(
        LinearViscosity(1.0e3),
        ParallelModel(PowerLawViscosity(4.0e-4, 3), LinearViscosity(1.0)),
    )
    εᵢⱼ = (1.0e-3, -1.0e-3, 0.0)
    εII = second_invariant_2D(εᵢⱼ)

    ks = x_keys(c)
    @test ks === (:τ, :ε)

    g = initial_guess_x(c, (; ε = εII, θ = 0.0), (; τ = 0.0), (;))
    @test all(isfinite, g)
    @test !iszero(g[2])
    @test g[2] == εII

    x = solve(c, g, (; ε = εᵢⱼ, θ = 0.0), (;); xnorm0 = normalisation_x(c, 1.0e6, εII))
    @test all(isfinite, x)
    @test x ≈ [0.011254510220682246, 0.000994372744889659] rtol = 1.0e-8

    # A tensor-valued vars.ε is reduced to its invariant for the seed.
    gt = initial_guess_x(c, (; ε = εᵢⱼ, θ = 0.0), (; τ = 0.0), (;))
    @test gt[2] == εII

    # A nonzero stress guess still produces the stress-based estimate; the
    # strain-rate fallback applies only where that estimate collapses.
    gτ = initial_guess_x(c, (; ε = εII, θ = 0.0), (; τ = 1.0e-2), (;))
    @test gτ[2] != εII
    @test gτ[2] > 0
end

@testset "solve reports non-convergence" begin
    c = SeriesModel(
        LinearViscosity(1.0e3),
        ParallelModel(PowerLawViscosity(4.0e-4, 3), LinearViscosity(1.0)),
    )
    εᵢⱼ = (1.0e-3, -1.0e-3, 0.0)
    εII = second_invariant_2D(εᵢⱼ)
    vars = (; ε = εᵢⱼ, θ = 0.0)
    xnorm = normalisation_x(c, 1.0e6, εII)

    # A zero branch strain rate is the singular seed: the residual becomes NaN,
    # which must surface as an error rather than as a NaN solution vector.
    @test_throws "did not converge" solve(c, SA[2.0, 0.0], vars, (;); xnorm0 = xnorm)
    err = try
        solve(c, SA[2.0, 0.0], vars, (;); xnorm0 = xnorm)
    catch e
        e
    end
    @test err isa NonConvergenceError
    @test isnan(err.residual)
    @test occursin("effective viscosity diverges", sprint(showerror, err))

    # Stopping short of convergence is also an error, with a different message.
    # Unreachable tolerances force the iteration to run out at `itermax`.
    err2 = try
        solve(c, initial_guess_x(c, (; ε = εII, θ = 0.0), (; τ = 0.0), (;)), vars, (;);
            xnorm0 = xnorm, atol = 0.0, rtol = 0.0, itermax = 5)
    catch e
        e
    end
    @test err2 isa NonConvergenceError
    @test isfinite(err2.residual)
    @test occursin("iteration limit", sprint(showerror, err2))
end
