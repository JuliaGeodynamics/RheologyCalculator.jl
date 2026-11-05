using RheologyCalculator, Test, StaticArrays, ForwardDiff
using RheologyCalculator.RheologyModels
using RheologyCalculator.RheologyModels: ModCamClay
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

@testset "initial guess with a zero-viscosity branch element" begin
    # compute_strain_rate of LinearViscosity(0) at τ = 0 is 0/0, which used to
    # propagate into the branch strain-rate seed.
    c = SeriesModel(
        IncompressibleElasticity(5.0),
        ParallelModel(PowerLawViscosity(4.0e-4, 3), LinearViscosity(0.0)),
    )
    εII = 1.0e-2
    g = initial_guess_x(c, (; ε = εII), (; τ = 0.0), (; dt = 0.1, τ0 = (0.0,)))
    @test all(isfinite, g)
    @test g[2] == εII
end

@testset "initial guess differentiated at zero strain rate" begin
    # A Dual strain rate with zero primal and nonzero partials must take the
    # zero branch of safe_inv, otherwise the seed carries Inf/NaN partials.
    @test RheologyCalculator.safe_inv(ForwardDiff.Dual(0.0, 4.0)) == ForwardDiff.Dual(0.0, 0.0)
    @test RheologyCalculator.safe_inv_one(ForwardDiff.Dual(0.0, 4.0)) == ForwardDiff.Dual(1.0, 0.0)

    c = SeriesModel(LinearViscosity(2.0), IncompressibleElasticity(5.0))
    others = (; dt = 0.1, τ0 = (0.7,))
    x0(ε) = initial_guess_x(c, (; ε = ε), (; τ = 0.0), others)
    @test isfinite(ForwardDiff.derivative(ε -> x0(ε)[1], 0.0))

    # dτ/dε through solve equals twice the series effective viscosity.
    τ(ε) = solve(c, x0(ε), (; ε = ε), others)[1]
    η_eff = inv(inv(2.0) + inv(5.0 * 0.1))
    @test ForwardDiff.derivative(τ, 0.0) ≈ 2 * η_eff
    @test ForwardDiff.derivative(τ, 1.0e-3) ≈ 2 * η_eff
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
        solve(
            c, initial_guess_x(c, (; ε = εII, θ = 0.0), (; τ = 0.0), (;)), vars, (;);
            xnorm0 = xnorm, atol = 0.0, rtol = 0.0, itermax = 5
        )
    catch e
        e
    end
    @test err2 isa NonConvergenceError
    @test isfinite(err2.residual)
    @test occursin("iteration limit", sprint(showerror, err2))

    @testset "floating-point stagnation is reported early" begin
        c = SeriesModel(
            LinearViscosity(1.0e23),
            Elasticity(1.0e10, 2.0e11),
            ModCamClay(; M = 0.9, N = 0.5, r = 1.0e8, β = 0.1, Pt = -1.0e5, η_vp = 1.0e20),
        )
        vars = (; ε = (0.0, -0.0, 0.0), θ = -7.0e-15)
        others = (; dt = 1.0e8, τ0 = ((0.0, 0.0, 0.0),), P0 = (1.3148e8,))
        x0 = SA[0.0, 0.0, 1.3148e8]
        xnorm = SA[7.0e-15, 1.0e10, 7.0e-15]

        err3 = try
            solve(c, x0, vars, others; xnorm0 = xnorm)
        catch e
            e
        end
        @test err3 isa NonConvergenceError
        @test err3.reason == :stagnation
        @test err3.iterations < 100
    end

    @testset "an iterate that hops by an ulp counts as stagnant" begin
        # Warm-started from its own solution without `xnorm0`, this model sits
        # at a residual of 7.3e-12, just above `atol`, and each Newton step moves
        # the iterate back and forth by an ulp or two.
        v1, v2 = LinearViscosity(5.0e19), LinearViscosity(1.0e20)
        c = SeriesModel(v1, ParallelModel(SeriesModel(v1, v2), v2), ParallelModel(v1, v2))
        vars = (; ε = 1.0e-15)
        sol = solve(c, initial_guess_x(c, vars, (; τ = 0.0), (;)), vars, (;))

        err4 = try
            solve(c, sol.x, vars, (;))
        catch e
            e
        end
        @test err4 isa NonConvergenceError
        @test err4.reason == :stagnation
        @test err4.iterations < 10
        @test err4.x ≈ sol.x rtol = 1.0e-14
    end

    @test RheologyCalculator.within_ulps(SA[1.0, 0.0], SA[nextfloat(1.0, 4), 0.0])
    @test !RheologyCalculator.within_ulps(SA[1.0, 0.0], SA[nextfloat(1.0, 5), 0.0])
    @test !RheologyCalculator.within_ulps(SA[1.0, 0.0], SA[1.0, 1.0e-300])
end
