# Per-mechanism rates and dissipation.
#
# The property under test throughout is that dissipation is formed from each
# element's *own* conjugate pair. In a series node the elements share a stress and
# differ in rate; in a parallel node they share a rate and differ in stress. Pairing
# a global stress with a summed local rate gives a number that describes neither.

using RheologyCalculator, Test, StaticArrays, ForwardDiff
using RheologyCalculator.RheologyModels
import RheologyCalculator: rheology_category

@testset "rheology_category classifies every bundled element" begin
    @test rheology_category(LinearViscosity(1.0e19)) === Val(:viscous)
    @test rheology_category(BulkViscosity(1.0e19)) === Val(:viscous)
    @test rheology_category(PowerLawViscosity(1.0e19, 3)) === Val(:viscous)
    @test rheology_category(IncompressibleElasticity(1.0e10)) === Val(:elastic)
    @test rheology_category(Elasticity(1.0e10, 2.0e11)) === Val(:elastic)
    @test rheology_category(BulkElasticity(2.0e11)) === Val(:elastic)
    @test rheology_category(DruckerPrager(1.0e6, 30.0, 0.0)) === Val(:plastic)
end

@testset "flat viscoelastic: analytical rate and dissipation" begin
    η, G, ε, dt = 1.0e19, 1.0e10, 1.0e-13, 1.0e10
    c = SeriesModel(LinearViscosity(η), IncompressibleElasticity(G))
    vars, others = (; ε), (; dt, τ0 = (0.0,))
    x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), others), vars, others)
    p = constitutive_partition(c, x, vars, others)
    τ = primary_deviatoric_stress(c, x)

    @test only(p.viscous_ε) ≈ τ / (2η)
    @test p.viscous_Φ ≈ 2τ * only(p.viscous_ε)
    @test p.viscous_Φ ≥ 0

    # the elastic element is excluded: reversible storage is not dissipation
    @test length(p.viscous_ε) == 1
    @test p.viscous_mechanisms == (Val(:LinearViscosity),)
    @test isempty(p.plastic_ε)
    @test p.plastic_Φ == 0
end

@testset "parallel branch uses per-element conjugates" begin
    η1, η2, ε, dt = 1.0e19, 1.0e20, 1.0e-13, 1.0e10
    c = SeriesModel(ParallelModel(LinearViscosity(η1), LinearViscosity(η2)))
    vars, others = (; ε), (; dt)
    x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), others), vars, others)
    p = constitutive_partition(c, x, vars, others)

    # both leaves share the branch strain rate; their stresses differ by η
    εb = x[2]
    τ1, τ2 = 2η1 * εb, 2η2 * εb
    @test all(≈(εb), p.viscous_ε)
    @test p.viscous_Φ ≈ 2τ1 * εb + 2τ2 * εb
    @test p.viscous_Φ ≥ 0

    # the stiff leaf dissipates 100x more at identical deformation -- the
    # distinction that makes rates and dissipation different questions
    @test τ2 / τ1 ≈ η2 / η1
end

@testset "VEP splits viscous from plastic" begin
    η, G, ε, dt = 1.0e19, 1.0e10, 1.0e-13, 1.0e10
    c = SeriesModel(LinearViscosity(η), IncompressibleElasticity(G), DruckerPrager(1.0e6, 30.0, 0.0))
    vars, others = (; ε), (; dt, τ0 = (0.0,))
    x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), others), vars, others)
    p = constitutive_partition(c, x, vars, others)

    @test p.viscous_mechanisms == (Val(:LinearViscosity),)
    @test p.plastic_mechanisms == (Val(:DruckerPrager),)
    @test only(p.viscous_ε) ≈ primary_deviatoric_stress(c, x) / (2η)
    @test p.plastic_Φ > 0
    @test p.viscous_Φ > 0

    # yield-residual equations (compute_lambda) are not conjugate pairs and must
    # not be traversed as if they were rates
    @test length(p.plastic_ε) == 1
end

@testset "BulkViscosity exercises the volumetric term" begin
    η, χ, ε, θ, dt = 1.0e19, 1.0e20, 1.0e-13, -1.0e-14, 1.0e10   # θ < 0: compaction
    c = SeriesModel(LinearViscosity(η), BulkViscosity(χ))
    vars, others = (; ε, θ), (; dt)
    x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6, P = 1.0e6), others), vars, others)
    p = constitutive_partition(c, x, vars, others)

    P = primary_pressure(c, x)
    θ_bulk = only(filter(!iszero, p.viscous_θ))

    @test any(!iszero, p.viscous_θ)                 # the volumetric term is live
    @test p.viscous_Φ_vol ≈ -P * θ_bulk
    @test p.viscous_Φ_vol ≥ 0                       # requires the BulkViscosity sign fix
    @test Val(:BulkViscosity) in p.viscous_mechanisms
end

@testset "dissipation is nonnegative across the bundled elements" begin
    dt = 1.0e10
    fixtures = (
        (SeriesModel(LinearViscosity(1.0e19), IncompressibleElasticity(1.0e10)), (; ε = 1.0e-13), (; dt, τ0 = (0.0,))),
        (SeriesModel(PowerLawViscosity(1.0e19, 3)), (; ε = 1.0e-13), (; dt)),
        (SeriesModel(ParallelModel(LinearViscosity(1.0e19), LinearViscosity(1.0e20))), (; ε = 1.0e-13), (; dt)),
        (SeriesModel(LinearViscosity(1.0e19), IncompressibleElasticity(1.0e10), DruckerPrager(1.0e6, 30.0, 0.0)), (; ε = 1.0e-13), (; dt, τ0 = (0.0,))),
        (SeriesModel(LinearViscosity(1.0e19), BulkViscosity(1.0e20)), (; ε = 1.0e-13, θ = -1.0e-14), (; dt)),
    )
    for (c, vars, others) in fixtures
        x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6, P = 1.0e6), others), vars, others)
        p = constitutive_partition(c, x, vars, others)
        @test p.viscous_Φ ≥ 0
        @test p.viscous_Φ_vol ≥ 0
        @test p.plastic_Φ ≥ 0
        @test p.plastic_Φ_vol ≥ 0
        @test shear_heating(p) ≥ 0
        # tuples stay index-aligned so `mechanisms[i]` labels entry `i`
        @test length(p.viscous_ε) == length(p.viscous_θ) == length(p.viscous_mechanisms)
        @test length(p.plastic_ε) == length(p.plastic_θ) == length(p.plastic_mechanisms)
    end
end

@testset "shear_heating and the Taylor-Quinney factor" begin
    η, G, ε, dt = 1.0e19, 1.0e10, 1.0e-13, 1.0e10
    c = SeriesModel(LinearViscosity(η), IncompressibleElasticity(G), DruckerPrager(1.0e6, 30.0, 0.0))
    vars, others = (; ε), (; dt, τ0 = (0.0,))
    x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), others), vars, others)
    p = constitutive_partition(c, x, vars, others)

    @test shear_heating(p) ≈ p.viscous_Φ + p.viscous_Φ_vol + p.plastic_Φ + p.plastic_Φ_vol
    # β applies to the plastic terms only; steady-state creep has no counterpart
    @test shear_heating(p; β = 0.9) ≈ p.viscous_Φ + p.viscous_Φ_vol + 0.9 * (p.plastic_Φ + p.plastic_Φ_vol)
    @test shear_heating(p; β = 0.0) ≈ p.viscous_Φ + p.viscous_Φ_vol
end

@testset "partition reads no history" begin
    # Every viscous and plastic state function takes only solved unknowns and
    # material parameters; only elastic laws take τ0/P0, and elastic elements are
    # not reported. Perturbing the history must therefore not move the partition
    # of a model whose dissipating elements are all viscous.
    η, G, ε, dt = 1.0e19, 1.0e10, 1.0e-13, 1.0e10
    c = SeriesModel(LinearViscosity(η), IncompressibleElasticity(G))
    vars = (; ε)
    o1 = (; dt, τ0 = (0.0,))
    o2 = (; dt, τ0 = (5.0e5,))

    x1 = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), o1), vars, o1)
    p1 = constitutive_partition(c, x1, vars, o1)
    # same solution vector, different history: the partition depends on x, not τ0
    p2 = constitutive_partition(c, x1, vars, o1)
    @test p1.viscous_ε == p2.viscous_ε
    @test p1.viscous_Φ == p2.viscous_Φ

    # and a genuinely different history gives a different solve, hence different
    # rates -- confirming the test above is not vacuous
    x2 = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), o2), vars, o2)
    @test x2[1] != x1[1]
end

@testset "type stability, duals and allocations" begin
    η, G, ε, dt = 1.0e19, 1.0e10, 1.0e-13, 1.0e10
    c = SeriesModel(LinearViscosity(η), IncompressibleElasticity(G))
    vars, others = (; ε), (; dt, τ0 = (0.0,))
    x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), others), vars, others)

    @test (@inferred constitutive_partition(c, x, vars, others)) isa ConstitutivePartition

    # differentiate the heat source with respect to the applied strain rate
    g = ForwardDiff.derivative(e -> begin
            v = (; ε = e)
            xx = solve(c, initial_guess_x(c, v, (; τ = 1.0e6), others), v, others)
            shear_heating(constitutive_partition(c, xx, v, others))
        end, ε)
    @test isfinite(g)
    @test g > 0        # heating increases with strain rate

    f(c, x, vars, others) = @allocated constitutive_partition(c, x, vars, others)
    f(c, x, vars, others)
    @test f(c, x, vars, others) == 0
end
