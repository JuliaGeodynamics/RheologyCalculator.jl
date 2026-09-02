using Test
using StaticArrays
import RheologyCalculator: AbstractRheology, compute_viscosity,
    compute_viscosity_series, compute_viscosity_parallel, correct_xnorm,
    _checked_η_KV, effective_strain_rate_correction, extract_local_kwargs

struct SimplificationTestRheology <: AbstractRheology end
compute_viscosity(::SimplificationTestRheology; kwargs...) = 3.0

# An element with no viscosity method at all: `compute_viscosity` falls back to
# 0.0, so a branch containing only this and a spring has no stiffness once `dt`
# is absent.
struct NoViscosityRheology <: AbstractRheology end

@testset "simplification-plan regressions" begin
    @test compute_viscosity_series(SimplificationTestRheology()) == 3.0
    @test compute_viscosity_parallel(SimplificationTestRheology()) == 3.0

    ltp = RheologyCalculator.RheologyModels.LTPViscosity(6.2e-13, 76.0, 1.8e9, 3.4e9)
    @test !isnothing(compute_viscosity_parallel(ltp))

    x = @SVector [1.0, 2.0]
    @test correct_xnorm(x, @SVector [3.0, 4.0]) == @SVector [3.0, 4.0]
    @test correct_xnorm(x, nothing) == @SVector [1.0, 1.0]
    @test_throws DimensionMismatch correct_xnorm(x, @SVector [1.0])
    @test_throws MethodError correct_xnorm(x, "invalid")

    # A branch that carries elastic history cannot have a zero effective
    # Kelvin-Voigt viscosity; `Inf` is legitimate and must still pass.
    elastic = IncompressibleElasticity(1.0e10)
    viscous = LinearViscosity(1.0e20)
    plastic = DruckerPrager(1.0e6, 30.0, 0.0)
    kv_args = (; ε = 1.0e-14, dt = 1.0e10, P = 1.0e6)

    @test _checked_η_KV((elastic, viscous), (), kv_args) ≈ 1.0e10 * 1.0e10 + 1.0e20
    @test _checked_η_KV((elastic, plastic), (), kv_args) == Inf
    @test_throws ArgumentError _checked_η_KV((elastic,), (), (; ε = 1.0e-14, dt = 0.0))
    @test_throws "must supply a nonzero `dt`" _checked_η_KV((elastic,), (), (; ε = 1.0e-14, dt = 0.0))

    # The same guard reached through the public correction entry point: `others`
    # without `dt` makes the spring's G*dt vanish, and the branch's only other
    # element contributes nothing.
    c_no_stiffness = SeriesModel(LinearViscosity(1.0e22), ParallelModel(elastic, NoViscosityRheology()))
    @test_throws ArgumentError effective_strain_rate_correction(
        c_no_stiffness, 1.0e-14, (2.0e6,), (; P0 = (0.0,))
    )
    @test effective_strain_rate_correction(
        c_no_stiffness, 1.0e-14, (2.0e6,), (; dt = 1.0e10, P0 = (0.0,))
    ) ≈ 2.0e6 / (2 * 1.0e10 * 1.0e10)

    # A history field must carry one entry per element claiming it; an element
    # index past the end is a malformed `others`, not a value to be guessed.
    others_hist = (; dt = 1.0e10, τ0 = (1.1, 3.0), d = (4, 2))
    @test extract_local_kwargs(others_hist, (:τ0,), 2) == (; dt = 1.0e10, τ0 = 3.0, d = 4)
    @test extract_local_kwargs(others_hist, (:d,), 2) == (; dt = 1.0e10, τ0 = 1.1, d = 2)
    @test_throws BoundsError extract_local_kwargs((; τ0 = (1.1,)), (:τ0,), 2)

    # Two springs in the same composite each claim τ0 and P0, so both tuples
    # need two entries.
    c_two_springs = SeriesModel(elastic, ParallelModel(elastic, viscous))
    @test_throws BoundsError initial_guess_x(
        c_two_springs, (; ε = 1.0e-15), (; τ = 1.0e2), (; dt = 1.0e10, τ0 = (0.0, 0.0), P0 = (0.0,))
    )

    # Element files under src/rheology/ extend the core generics; they do not
    # define same-named functions local to RheologyModels. A local definition
    # would leave the core generic on its fallback for that element, which is
    # how a missing compute_viscosity_series silently halves a Kelvin-Voigt
    # aggregate.
    @testset "element files extend the core generics" begin
        extended = (
            :series_state_functions, :parallel_state_functions, :_isvolumetric, :isvolumetric,
            :compute_strain_rate, :compute_stress, :compute_pressure,
            :compute_volumetric_strain_rate, :compute_plastic_strain_rate,
            :compute_plastic_stress, :compute_volumetric_plastic_strain_rate,
            :compute_lambda, :compute_lambda_parallel,
            :compute_viscosity, :compute_viscosity_series, :compute_viscosity_parallel,
        )
        for name in extended
            @test getfield(RheologyCalculator.RheologyModels, name) ===
                getfield(RheologyCalculator, name)
        end
    end
end
