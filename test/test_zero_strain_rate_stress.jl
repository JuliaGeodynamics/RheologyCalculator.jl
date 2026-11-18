# A zero strain rate with no stress history must solve to exactly τ = 0.
using Test, ForwardDiff
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "zero strain rate gives τ = 0" begin
    zero_ε(::Val{:scalar}) = 0.0
    zero_ε(::Val{:tensor}) = (0.0, 0.0, 0.0)
    zero_τ0(::Val{:scalar}) = (0.0,)
    zero_τ0(::Val{:tensor}) = ((0.0, 0.0, 0.0),)

    models = (
        "Maxwell" => (SeriesModel(LinearViscosity(1.0e21), IncompressibleElasticity(1.0e10)), (; dt = 1.0e9), (; τ = 0.0)),
        "power law + elasticity" => (SeriesModel(PowerLawViscosity(1.0e19, 3), IncompressibleElasticity(1.0e10)), (; dt = 1.0e9), (; τ = 0.0)),
        "VEP Drucker-Prager" => (
            SeriesModel(LinearViscosity(1.0e21), IncompressibleElasticity(1.0e10), DruckerPrager(10.0e6, 30, 0)),
            (; dt = 1.0e9, P = 1.0e6, P0 = (0.0,)), (; τ = 0.0, λ = 0.0),
        ),
    )

    for (name, (c, extra, args)) in models, form in (Val(:scalar), Val(:tensor))
        @testset "$name, $(typeof(form).parameters[1])" begin
            vars = (; ε = zero_ε(form))
            others = merge(extra, (; τ0 = zero_τ0(form)))
            x = initial_guess_x(c, vars, args, others)
            sol = solve(c, x, vars, others)
            @test sol.x[1] === 0.0
        end
    end

    @testset "values away from zero are unchanged" begin
        c = SeriesModel(LinearViscosity(1.0), IncompressibleElasticity(10.0))
        others = (; dt = 0.1, τ0 = (0.0,))
        τ = solve(c, initial_guess_x(c, (; ε = 2.0), (; τ = 0.0), others), (; ε = 2.0), others).x[1]
        # Maxwell: ε = τ/(2η) + τ/(2G dt)
        @test τ ≈ 2.0 / (1 / 2 + 1 / (2 * 10.0 * 0.1)) rtol = 1.0e-10
    end

    @testset "derivative through solve at ε = 0 is finite" begin
        c = SeriesModel(LinearViscosity(1.0e21), IncompressibleElasticity(1.0e10))
        others = (; dt = 1.0e9, τ0 = (0.0,))
        τ(ε) = solve(c, initial_guess_x(c, (; ε), (; τ = 0.0), others), (; ε), others).x[1]
        d0 = ForwardDiff.derivative(τ, 0.0)
        @test isfinite(d0)
        # dτ/dε = 1 / (1/(2η) + 1/(2G dt))
        @test d0 ≈ 1 / (1 / (2 * 1.0e21) + 1 / (2 * 1.0e10 * 1.0e9)) rtol = 1.0e-8
    end
end
