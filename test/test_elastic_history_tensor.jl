include(joinpath(@__DIR__, "reference", "tensor_reference.jl"))

using .TensorReference: ReferenceModel, reference_step, initial_guess, tensor_invariant
using RheologyCalculator.RheologyModels: elastic_stress_history_2D, elastic_stress_history_3D

# Time loops of `solve` + `elastic_stress_history_*` compared step by step with
# the full-tensor reference solver in `reference/tensor_reference.jl`.

module ElasticHistoryCases

using RheologyCalculator, RheologyCalculator.RheologyModels

const dt = 1.0e10
const nsteps = 12
const εII = 1.0e-14
# Characteristic scales for the reference solver's unknowns.
const τc = 1.0e7
const εc = εII

# Relaxation times η/G ~ 5e10 s, comparable to a few steps.
const η₁ = LinearViscosity(1.0e21)
const η₂ = LinearViscosity(5.0e20)
const η₃ = LinearViscosity(8.0e20)
const η₄ = LinearViscosity(3.0e20)
const G₁ = IncompressibleElasticity(1.0e10)
const G₂ = IncompressibleElasticity(2.5e10)
const Ga = IncompressibleElasticity(1.0e10)
const Gb = IncompressibleElasticity(4.0e10)
# ε = τⁿ/(2η): effective viscosity ~5e20 Pa s at τ ~ 1e7 Pa.
const PL = PowerLawViscosity(2.5e34, 3)

const linear_topologies = (
    Maxwell = SeriesModel(η₁, G₁),
    KV = SeriesModel(η₁, ParallelModel(η₂, G₁)),
    mixed = SeriesModel(η₁, ParallelModel(η₂, SeriesModel(η₃, G₁))),
    Burgers = SeriesModel(η₁, G₁, ParallelModel(η₂, G₂)),
    two_blocks = SeriesModel(η₁, ParallelModel(η₂, G₁), ParallelModel(η₃, SeriesModel(η₄, G₂))),
    spring_and_subbranch = SeriesModel(η₁, ParallelModel(G₁, SeriesModel(η₃, G₂))),
    two_springs_subbranch = SeriesModel(η₁, ParallelModel(η₂, SeriesModel(η₃, Ga, Gb))),
)

const nonlinear_topologies = (
    PL_direct = SeriesModel(η₁, ParallelModel(PL, G₁)),
    PL_subbranch = SeriesModel(η₁, ParallelModel(η₂, SeriesModel(PL, G₁))),
    PL_outer = SeriesModel(PL, ParallelModel(η₂, G₁)),
)

const pure_shear_2D = (1.0, -1.0, 0.0)
const simple_shear_2D = (0.0, 0.0, 1.0)
const dir_a_3D = (1.0, -0.4, -0.6, 0.3, 0.2, 0.5)
const dir_b_3D = (-0.2, 0.7, -0.5, -0.4, 0.6, 0.1)

_tensor(d) = map(di -> εII * di / RheologyCalculator.second_invariant(d...), d)
_flip(t) = map(-, t)

# Each load case maps step `n` (1-based) to the applied strain rate.
const half = nsteps ÷ 2
const load_cases = (
    pure_shear = n -> _tensor(pure_shear_2D),
    reversal = n -> n ≤ half ? _tensor(pure_shear_2D) : _flip(_tensor(pure_shear_2D)),
    rotation = n -> n ≤ half ? _tensor(pure_shear_2D) : _tensor(simple_shear_2D),
    three_D = n -> n ≤ half ? _tensor(dir_a_3D) : _tensor(dir_b_3D),
    scalar = n -> εII,
    scalar_reversal = n -> n ≤ half ? εII : -εII,
)

end # module

"""
    compare_histories(c, load; feed)

Run `nsteps` steps of the package and of the reference for composite `c` under
the strain-rate schedule `load`. Returns the largest relative errors over all
steps: `τII` (outer stress invariant, signed for scalar ε) and `springs` (every
spring tensor from `elastic_stress_history_*`).

With `feed = :package` each side advances with its own spring history, which is
the end-to-end time loop. With `feed = :reference` the package solves every step
from the reference history, so `τII` measures the solve alone and `springs` the
history recovery alone.
"""
function compare_histories(c, load; feed::Symbol)
    feed in (:package, :reference) || throw(ArgumentError("feed must be :package or :reference"))
    (; dt, nsteps, τc, εc) = ElasticHistoryCases
    m = ReferenceModel(c)
    ε1 = load(1)
    nc = ε1 isa Number ? 1 : length(ε1)
    as_tuple(t) = t isa Number ? (t,) : t
    zero_t = ntuple(_ -> 0.0, nc)
    nsprings = TensorReference.n_springs(m)
    hist_pkg = ntuple(_ -> ε1 isa Number ? 0.0 : zero_t, nsprings)
    hist_ref = ntuple(_ -> zero_t, nsprings)
    history = nc == 6 ? elastic_stress_history_3D : elastic_stress_history_2D
    y = initial_guess(m, as_tuple(ε1))
    x = nothing
    xnorm0 = normalisation_x(c, τc, εc)
    err_τ = 0.0
    err_springs = 0.0
    for n in 1:nsteps
        ε = load(n)
        τ_ref, springs_ref, y = reference_step(m, as_tuple(ε), hist_ref, dt; guess = y, τc, εc)

        τ0 = feed === :package ? hist_pkg : hist_ref
        others = (; dt, τ0)
        vars = (; ε, θ = 0.0)
        x0 = x === nothing ? initial_guess_x(c, vars, (; τ = τc, P = 0.0), others) : x
        x = solve(c, x0, vars, others; xnorm0, itermax = 200).x
        τII_ref = nc == 1 ? τ_ref[1] : tensor_invariant(τ_ref)
        err_τ = max(err_τ, abs(x[1] - τII_ref) / abs(τII_ref))

        springs_pkg = history(c, x, ε, τ0, others)
        for (sp, sr) in zip(springs_pkg, springs_ref)
            scale = max(tensor_invariant(sr), 1.0e-6 * τc)
            err_springs = max(err_springs, tensor_invariant(map(-, as_tuple(sp), sr)) / scale)
        end
        hist_pkg = springs_pkg
        hist_ref = springs_ref
    end
    return (; τII = err_τ, springs = err_springs)
end

const history_rtol = 1.0e-9

@testset "full-tensor reference: Maxwell closed form" begin
    # Backward Euler for S(η, G): τ = η*(τᵒ + 2GΔt ε), η* = η/(η + GΔt).
    (; dt, τc, εc, η₁, G₁) = ElasticHistoryCases
    m = ReferenceModel(SeriesModel(η₁, G₁))
    ηstar = η₁.η / (η₁.η + G₁.G * dt)
    τᵒ = (3.0e6, -1.0e6, 2.0e6)
    ε = (0.0, 0.0, 1.0e-14)
    τ, springs, _ = reference_step(m, ε, (τᵒ,), dt; guess = initial_guess(m, ε), τc, εc)
    expected = map((a, b) -> ηstar * (a + 2 * G₁.G * dt * b), τᵒ, ε)
    @test all(isapprox.(τ, expected; rtol = 1.0e-12))
    @test springs == (τ,)
end

@testset "zero assembled strain rate ε + H" begin
    # The invariant of ε + H has no derivative at zero; the solve must still
    # return a finite solution there, both from rest and when the history
    # cancels the applied strain rate.
    (; dt, τc, εc, η₁, η₂, G₁) = ElasticHistoryCases
    ε = (1.0e-14, -1.0e-14, 0.0)
    zero_t = (0.0, 0.0, 0.0)
    for (c, cancelling_τ0) in (
            (SeriesModel(η₁, G₁), map(e -> -2 * G₁.G * dt * e, ε)),
            (SeriesModel(η₁, ParallelModel(η₂, G₁)), map(e -> -2 * (η₂.η + G₁.G * dt) * e, ε)),
        )
        xnorm0 = normalisation_x(c, τc, εc)
        for (vars, others) in (
                ((; ε = zero_t, θ = 0.0), (; dt, τ0 = (zero_t,))),
                ((; ε, θ = 0.0), (; dt, τ0 = (cancelling_τ0,))),
            )
            x0 = initial_guess_x(c, vars, (; τ = τc, P = 0.0), others)
            x = solve(c, x0, vars, others; xnorm0).x
            @test all(isfinite, x)
            @test abs(x[1]) < 1.0e-9 * τc
        end
    end
end

@testset "elastic history vs full-tensor reference" begin
    C = ElasticHistoryCases
    for (tn, c) in pairs(merge(C.linear_topologies, C.nonlinear_topologies)), (ln, load) in pairs(C.load_cases)
        @testset "$tn / $ln" begin
            fed = compare_histories(c, load; feed = :reference)
            own = compare_histories(c, load; feed = :package)
            @test fed.τII < history_rtol
            @test fed.springs < history_rtol
            @test max(own.τII, own.springs) < history_rtol
        end
    end
end

@testset "closure unknowns of blocks carrying elastic history" begin
    C = ElasticHistoryCases
    # Each closure unknown is a local equation of its block, ahead of the
    # sub-branches; blocks of linear elements keep their layout.
    @test x_keys(C.linear_topologies.KV) == (:τ, :ε)
    @test x_keys(C.linear_topologies.mixed) == (:τ, :ε, :τ)
    @test x_keys(C.nonlinear_topologies.PL_direct) == (:τ, :ε, :εII_block)
    @test x_keys(C.nonlinear_topologies.PL_subbranch) == (:τ, :ε, :τII_subbranch, :τ)
    @test x_keys(C.nonlinear_topologies.PL_outer) == (:τ, :ε)
    c_both = SeriesModel(C.η₁, ParallelModel(C.PL, SeriesModel(C.PL, C.G₁), SeriesModel(C.η₃, C.G₂)))
    @test x_keys(c_both) == (:τ, :ε, :εII_block, :τII_subbranch, :τ, :τ)

    # Closure unknowns are invariants, bounded for a tensor strain rate and
    # signed for a scalar one. Their residuals are scaled like the block
    # (strain rate) and sub-branch (stress) residuals.
    c = C.nonlinear_topologies.PL_subbranch
    @test RheologyCalculator.branch_strain_rate_mask(c, (1.0, -1.0, 0.0)) == [false, true, true, false]
    @test !any(RheologyCalculator.branch_strain_rate_mask(c, 1.0))
    @test normalisation_x(c, 2.0, 3.0) == [3.0, 2.0, 2.0, 3.0]

    # At the solution each closure unknown equals the invariant of the physical
    # tensor of the reference solve.
    (; dt, τc, εc) = C
    m = ReferenceModel(c)
    ε = C._tensor(C.pure_shear_2D)
    τ0 = ((2.0e6, -1.0e6, 1.5e6),)
    _, springs_ref, _ = reference_step(m, ε, τ0, dt; guess = initial_guess(m, ε), τc, εc)
    others = (; dt, τ0)
    vars = (; ε, θ = 0.0)
    x0 = initial_guess_x(c, vars, (; τ = τc, P = 0.0), others)
    x = solve(c, x0, vars, others; xnorm0 = normalisation_x(c, τc, εc)).x
    # The sub-branch carries the stress of its only spring.
    @test x[3] ≈ tensor_invariant(springs_ref[1]) rtol = 1.0e-9
end

@testset "spring stress recovery" begin
    C = ElasticHistoryCases
    (; dt, τc, εc) = C
    c = C.linear_topologies.two_blocks
    ε = C._tensor(C.dir_a_3D)
    τ0 = ((2.0e6, -1.0e6, -1.0e6, 5.0e5, 0.0, 1.0e6), (-1.0e6, 3.0e6, -2.0e6, 0.0, 1.0e6, -5.0e5))
    others = (; dt, τ0)
    vars = (; ε, θ = 0.0)
    x0 = initial_guess_x(c, vars, (; τ = τc, P = 0.0), others)
    x = solve(c, x0, vars, others; xnorm0 = normalisation_x(c, τc, εc)).x
    springs = elastic_stress_history(c, x, ε, τ0, others)
    @test length(springs) == length(τ0)
    @test springs == elastic_stress_history_3D(c, x, ε, τ0, others)
    @test compute_stress_elastic(c, x, ε, others) == map(tensor_invariant, springs)

    # A tensor history makes the block springs non-coaxial with the outer
    # stress, which the methods without ε cannot recover.
    @test_throws "use `compute_stress_elastic(c, x, ε, others)`" compute_stress_elastic(c, x, others)
    @test_throws "pass the solver vector `x`" elastic_stress_history_3D(c, x[1], ε, τ0, others)

    # A scalar history has no direction: both methods give the signed stresses.
    τ0s = (2.0e6, -1.0e6)
    others_s = (; dt, τ0 = τ0s)
    vars_s = (; ε = -C.εII, θ = 0.0)
    xs = solve(c, x0, vars_s, others_s; xnorm0 = normalisation_x(c, τc, εc)).x
    @test compute_stress_elastic(c, xs, others_s) == compute_stress_elastic(c, xs, vars_s.ε, others_s)
end

@testset "component partition inside blocks carrying elastic history" begin
    # Elements of a block are reported at the physical block strain rate and
    # sub-branch stress, not at the shifted unknowns of the solver.
    C = ElasticHistoryCases
    (; dt, τc, εc, η₁, η₂, G₁, PL) = C
    ε = C._tensor(C.pure_shear_2D)
    τ0 = ((3.0e6, 1.0e6, -2.0e6),)
    others = (; dt, τ0)
    vars = (; ε, θ = 0.0)
    for (c, block_elements) in (
            (C.linear_topologies.KV, (η₂,)),
            (C.nonlinear_topologies.PL_direct, (PL,)),
            (C.nonlinear_topologies.PL_subbranch, (η₂, PL)),
        )
        m = ReferenceModel(c)
        τ_ref, springs_ref, _ = reference_step(m, ε, τ0, dt; guess = initial_guess(m, ε), τc, εc)
        x0 = initial_guess_x(c, vars, (; τ = τc, P = 0.0), others)
        x = solve(c, x0, vars, others; xnorm0 = normalisation_x(c, τc, εc)).x
        p = component_partition(c, x, vars, others)
        @test p.viscous_elements == (η₁, block_elements...)
        # outer dashpot: the outer stress
        @test p.viscous_τ[1] ≈ tensor_invariant(τ_ref) rtol = 1.0e-9
        if c === C.nonlinear_topologies.PL_subbranch
            # direct dashpot of the block: the block strain rate, from the
            # stress balance τ = 2η₂ ε̇_b + τ_sub
            ε_block = tensor_invariant(map((t, s) -> (t - s) / (2 * η₂.η), τ_ref, springs_ref[1]))
            @test p.viscous_ε[2] ≈ ε_block rtol = 1.0e-9
            # power law of the sub-branch: the stress of its spring
            τ_sub = tensor_invariant(springs_ref[1])
            @test p.viscous_τ[3] ≈ τ_sub rtol = 1.0e-9
            @test p.viscous_ε[3] ≈ RheologyCalculator.compute_strain_rate(PL; τ = τ_sub) rtol = 1.0e-9
        else
            # direct element of the block: the block strain rate, from its spring
            ε_block = tensor_invariant(map((s, s0) -> (s - s0) / (2 * G₁.G * dt), springs_ref[1], τ0[1]))
            @test p.viscous_ε[2] ≈ ε_block rtol = 1.0e-9
            @test p.viscous_τ[2] ≈ RheologyCalculator.compute_stress(block_elements[1]; ε = ε_block) rtol = 1.0e-9
        end
    end
end
