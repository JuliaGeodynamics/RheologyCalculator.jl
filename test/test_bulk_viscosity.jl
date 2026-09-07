# Volumetric sign convention.
#
# `P` is positive in compression, `θ` positive in dilation. Compressing a material
# must therefore compact it, and the two must carry opposite signs. `BulkViscosity`
# previously returned `+P/χ`, which asserted that compression causes dilation and
# made the volumetric dissipation `-P*θ = -P^2/χ` negative -- a viscous element
# generating energy. It had no test coverage: `test_jacobians.jl` bound a
# `BulkViscosity` and never used it, and a Jacobian check cannot catch a sign
# error anyway, since a sign-flipped law is still self-consistent under AD.

using RheologyCalculator, Test
using RheologyCalculator.RheologyModels
import RheologyCalculator: compute_pressure, compute_volumetric_strain_rate

@testset "BulkViscosity volumetric sign" begin
    χ = 1.0e19
    bv = BulkViscosity(χ)

    # compression compacts, dilation puts the material in tension
    @test compute_volumetric_strain_rate(bv; P = 1.0e6) < 0
    @test compute_volumetric_strain_rate(bv; P = -1.0e6) > 0
    @test compute_pressure(bv; θ = 1.0e-15) < 0
    @test compute_pressure(bv; θ = -1.0e-15) > 0

    # the two directions are consistent inverses
    P = 1.0e6
    θ = compute_volumetric_strain_rate(bv; P)
    @test compute_pressure(bv; θ) ≈ P

    # volumetric dissipation -P*θ must be non-negative for any pressure
    for P in (-1.0e7, -1.0e3, 0.0, 1.0e3, 1.0e7)
        θ = compute_volumetric_strain_rate(bv; P)
        @test -P * θ ≥ 0
    end
end

@testset "BulkViscosity agrees in sign with BulkElasticity" begin
    bv = BulkViscosity(1.0e19)
    be = BulkElasticity(2.0e11)
    dt = 1.0e10

    # Both are volumetric laws under the same convention, so a compressive
    # pressure must compact in both. This is the comparison that would have
    # caught the original inversion.
    P = 1.0e6
    @test sign(compute_volumetric_strain_rate(bv; P)) ==
        sign(compute_volumetric_strain_rate(be; P, P0 = 0.0, dt))
end

@testset "BulkViscosity solves inside a composite" begin
    ε, θ, dt = 1.0e-15, -1.0e-15, 1.0e10   # compaction
    c = SeriesModel(LinearViscosity(1.0e19), BulkViscosity(1.0e20))
    @test x_keys(c) == (:τ, :P)

    vars = (; ε, θ)
    others = (; dt)
    x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e3, P = 1.0e5), others), vars, others)

    # compaction (θ < 0) must produce a compressive (positive) pressure
    @test primary_pressure(c, x) > 0
    @test primary_deviatoric_stress(c, x) > 0
end
