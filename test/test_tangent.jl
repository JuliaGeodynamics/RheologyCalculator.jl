using RheologyCalculator, Test, StaticArrays
using RheologyCalculator.RheologyModels
import ForwardDiff

# Central difference of the whole converged solve, used as the reference the
# analytical tangent must reproduce.
function fd_tangent(c, vars, others, xnorm; δ = 1.0e-6)
    εII = vars.ε
    τ(e) = begin
        v = merge(vars, (; ε = e))
        x = solve(c, initial_guess_x(c, v, (; τ = 0.0), others), v, others; xnorm0 = xnorm)
        x[stress_index(c)]
    end
    return (τ(εII * (1 + δ)) - τ(εII * (1 - δ))) / (2 * εII * δ)
end

converged(c, vars, others, xnorm) =
    solve(c, initial_guess_x(c, vars, (; τ = 0.0), others), vars, others; xnorm0 = xnorm)

@testset "tangent" begin
    @testset "linear viscosity is its own tangent" begin
        # τII = 2ηεII, so dτII/dεII = 2η exactly.
        η = 1.0e3
        c = SeriesModel(LinearViscosity(η))
        vars = (; ε = 1.0e-3, θ = 0.0)
        xnorm = normalisation_x(c, 1.0e6, vars.ε)
        x = converged(c, vars, (;), xnorm)
        @test tangent(c, x, vars, (;)) ≈ 2η rtol = 1.0e-12
    end

    @testset "power law against the analytical derivative" begin
        # εII = τII^n/(2η) inverted gives τII = (2ηεII)^(1/n),
        # so dτII/dεII = (2η)^(1/n) εII^(1/n - 1) / n.
        η, n = 4.0e-4, 3
        c = SeriesModel(PowerLawViscosity(η, n))
        εII = 1.0e-3
        vars = (; ε = εII, θ = 0.0)
        xnorm = normalisation_x(c, 1.0e6, εII)
        x = converged(c, vars, (;), xnorm)
        expected = (2η)^(1 / n) * εII^(1 / n - 1) / n
        @test tangent(c, x, vars, (;)) ≈ expected rtol = 1.0e-10
    end

    @testset "nested parallel branch matches finite differences" begin
        # The composite from issue #41, whose reference value is 5.09918017.
        c = SeriesModel(
            ParallelModel(
                SeriesModel(PowerLawViscosity(4.0e-4, 3), LinearViscosity(1.0e3)),
                LinearViscosity(1.0),
            ),
        )
        εII = 1.0e-3
        vars = (; ε = εII, θ = 0.0)
        xnorm = normalisation_x(c, 1.0e6, εII)
        x = solve(c, SA[0.011, εII, 0.009], vars, (;); xnorm0 = xnorm)
        @test tangent(c, x, vars, (;)) ≈ 5.09918017 rtol = 1.0e-8
    end

    @testset "viscoelastic composite matches finite differences" begin
        # Exercises the elastic path: τ0 in `others`, and the pre-solve
        # correction of ε that `tangent` has to reproduce.
        c = SeriesModel(LinearViscosity(1.0e22), IncompressibleElasticity(1.0e10))
        εII = 1.0e-14
        vars = (; ε = εII, θ = 0.0)
        others = (; dt = 1.0e10, τ0 = (1.0e6,), P0 = (0.0,))
        xnorm = normalisation_x(c, 1.0e6, εII)
        x = converged(c, vars, others, xnorm)
        @test tangent(c, x, vars, others) ≈ fd_tangent(c, vars, others, xnorm) rtol = 1.0e-6
    end

    @testset "series viscous + power law matches finite differences" begin
        c = SeriesModel(LinearViscosity(1.0e3), PowerLawViscosity(4.0e-4, 3))
        εII = 1.0e-3
        vars = (; ε = εII, θ = 0.0)
        xnorm = normalisation_x(c, 1.0e6, εII)
        x = converged(c, vars, (;), xnorm)
        @test tangent(c, x, vars, (;)) ≈ fd_tangent(c, vars, (;), xnorm) rtol = 1.0e-6
    end

    @testset "residual differentiates with respect to vars.ε" begin
        # The residual is heterogeneous when only one equation carries a Dual;
        # subtract_elastic_correction must still accept it.
        c = SeriesModel(
            LinearViscosity(1.0e22),
            ParallelModel(LinearViscosity(1.0e21), IncompressibleElasticity(1.0e10)),
        )
        εII = 1.0e-14
        others = (; dt = 1.0e10, τ0 = (1.0e6,), P0 = (0.0,))
        vars = (; ε = εII, θ = 0.0)
        x = converged(c, vars, others, normalisation_x(c, 1.0e6, εII))
        dRdε = ForwardDiff.derivative(
            e -> compute_residual(c, x, merge(vars, (; ε = e)), others), εII
        )
        @test all(isfinite, dRdε)
        @test !iszero(dRdε)
    end

    @testset "stress_index" begin
        c = SeriesModel(LinearViscosity(1.0e22), Elasticity(1.0e10, 1.0e11))
        @test x_keys(c)[stress_index(c)] === :τ

        # A composite without a τ unknown reports so rather than returning a
        # silently wrong entry.
        @test_throws "no deviatoric stress unknown" RheologyCalculator._stress_index((:P, :θ))
    end
end
