# Semantic access to the local solution vector.
#
# The point of these queries is that `x_keys` labels an unknown by its *kind*, not
# its location, so a nested composite legitimately repeats a key. Selecting "the"
# stress by symbol is therefore wrong in general; it must be selected by topology.

using RheologyCalculator, Test, StaticArrays, ForwardDiff
using RheologyCalculator.RheologyModels

@testset "solution_values preserves every entry with a key" begin
    v, e = LinearViscosity(1.0e18), IncompressibleElasticity(1.0e10)
    # a Maxwell element inside a parallel branch => two distinct :τ unknowns
    c = SeriesModel(v, ParallelModel(SeriesModel(LinearViscosity(1.0e19), e), v))
    @test x_keys(c) == (:τ, :ε, :τ)

    x = SA[10.0, 20.0, 30.0]
    @test solution_values(c, x, Val(:τ)) == (10.0, 30.0)   # both, in order
    @test solution_values(c, x, Val(:ε)) == (20.0,)
    @test solution_values(c, x, Val(:P)) == ()             # absent key -> empty
end

@testset "primary_deviatoric_stress selects the global equation, not the first symbol" begin
    v, e = LinearViscosity(1.0e18), IncompressibleElasticity(1.0e10)
    ε, dt = 1.0e-14, 1.0e10
    η_soft, η_stiff = 1.0e19, 1.0e18

    # Two layouts, same physical primary stress. The nested one repeats :τ.
    flat = SeriesModel(v, e)
    nested = SeriesModel(v, ParallelModel(SeriesModel(LinearViscosity(η_soft), e), LinearViscosity(η_stiff)))

    for c in (flat, nested)
        others = (; dt, τ0 = (0.0,))
        sol = solve(c, initial_guess_x(c, (; ε), (; τ = 1.0), others), (; ε), others)
        # the global stress is the first entry for every topology the package builds
        @test primary_deviatoric_stress(c, sol) == sol[1]
    end

    # On the nested model the two :τ entries differ by orders of magnitude, so
    # picking the wrong one is not a rounding error.
    others = (; dt, τ0 = (0.0,))
    sol = solve(nested, initial_guess_x(nested, (; ε), (; τ = 1.0), others), (; ε), others)
    @test sol[1] != sol[3]
    @test primary_deviatoric_stress(nested, sol) == sol[1]

    # Works on a bare SVector too, not just whatever `solve` returns.
    @test primary_deviatoric_stress(nested, SA[1.0, 2.0, 3.0]) == 1.0
end

@testset "primary_pressure and its fallback" begin
    v = LinearViscosity(1.0e18)
    compressible = SeriesModel(v, Elasticity(1.0e10, 2.0e11))
    incompressible = SeriesModel(v, IncompressibleElasticity(1.0e10))

    @test x_keys(compressible) == (:τ, :P)
    @test primary_pressure(compressible, SA[10.0, 99.0]) == 99.0

    # No volumetric unknown: without a fallback this is `nothing`, with one it is
    # the caller's value. The fallback is policy, not solution data.
    @test x_keys(incompressible) == (:τ,)
    @test primary_pressure(incompressible, SA[10.0]) === nothing
    @test primary_pressure(incompressible, SA[10.0]; fallback = 1.0e6) == 1.0e6
end

@testset "plastic_multipliers" begin
    v, e = LinearViscosity(1.0e18), IncompressibleElasticity(1.0e10)
    plastic = SeriesModel(v, e, DruckerPrager(1.0e6, 30.0, 0.0))
    viscous = SeriesModel(v, e)

    @test x_keys(plastic) == (:τ, :λ)
    @test plastic_multipliers(plastic, SA[10.0, 0.5]) == (0.5,)
    @test plastic_multipliers(viscous, SA[10.0]) == ()      # empty, not an error
end

@testset "ambiguous by-name access is refused, semantic access is not" begin
    v, e = LinearViscosity(1.0e18), IncompressibleElasticity(1.0e10)
    c = SeriesModel(v, ParallelModel(SeriesModel(LinearViscosity(1.0e19), e), v))
    x = SA[10.0, 20.0, 30.0]
    # `solution_values` reports both; the semantic query still resolves uniquely.
    @test length(solution_values(c, x, Val(:τ))) == 2
    @test primary_deviatoric_stress(c, x) == 10.0
end

@testset "stress_index agrees with the topological query on every bundled layout" begin
    v = LinearViscosity(1.0e18)
    e = IncompressibleElasticity(1.0e10)
    models = (
        SeriesModel(v),
        SeriesModel(v, e),
        SeriesModel(v, Elasticity(1.0e10, 2.0e11)),
        SeriesModel(v, e, DruckerPrager(1.0e6, 30.0, 0.0)),
        SeriesModel(v, ParallelModel(LinearViscosity(1.0e19), e)),
        SeriesModel(v, ParallelModel(SeriesModel(LinearViscosity(1.0e19), e), v)),
        SeriesModel(ParallelModel(LinearViscosity(1.0e19), e), ParallelModel(LinearViscosity(1.0e20), e)),
    )
    for c in models
        # The replacement must not move any existing answer.
        @test stress_index(c) == primary_stress_index(c)
    end
end

@testset "queries are type stable, dual-compatible and allocation free" begin
    v, e = LinearViscosity(1.0e18), IncompressibleElasticity(1.0e10)
    c = SeriesModel(v, ParallelModel(SeriesModel(LinearViscosity(1.0e19), e), v))
    x = SA[10.0, 20.0, 30.0]

    @test (@inferred primary_deviatoric_stress(c, x)) == 10.0
    @test (@inferred solution_values(c, x, Val(:τ))) == (10.0, 30.0)
    @test (@inferred plastic_multipliers(c, x)) == ()

    # differentiating through the query must work
    g = ForwardDiff.derivative(t -> primary_deviatoric_stress(c, SA[t, 20.0, 30.0]), 10.0)
    @test g == 1.0

    f1(c, x) = @allocated primary_deviatoric_stress(c, x)
    f2(c, x) = @allocated solution_values(c, x, Val(:τ))
    f3(c, x) = @allocated plastic_multipliers(c, x)
    f1(c, x); f2(c, x); f3(c, x)          # warm
    @test f1(c, x) == 0
    @test f2(c, x) == 0
    @test f3(c, x) == 0
end
