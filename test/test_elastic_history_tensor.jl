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

An error from the package ends the comparison and reports `NaN` for what it
prevented: a failed `solve` (non-convergence, or a `DomainError` from an element
law evaluated outside its domain) makes both errors `NaN`; a failed history recovery
makes `springs` `NaN`, and also `τII` when `feed = :package`.
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
        x = try
            solve(c, x0, vars, others; xnorm0, itermax = 200).x
        catch err
            err isa Union{NonConvergenceError, DomainError} || rethrow()
            return (; τII = NaN, springs = NaN)
        end
        τII_ref = nc == 1 ? τ_ref[1] : tensor_invariant(τ_ref)
        err_τ = max(err_τ, abs(x[1] - τII_ref) / abs(τII_ref))

        springs_pkg = try
            history(c, x, ε, τ0, others)
        catch err
            err isa ArgumentError || rethrow()
            feed === :package && return (; τII = NaN, springs = NaN)
            err_springs = NaN
            nothing
        end
        if springs_pkg !== nothing
            for (sp, sr) in zip(springs_pkg, springs_ref)
                scale = max(tensor_invariant(sr), 1.0e-6 * τc)
                err_springs = max(err_springs, tensor_invariant(map(-, as_tuple(sp), sr)) / scale)
            end
            hist_pkg = springs_pkg
        end
        hist_ref = springs_ref
    end
    return (; τII = err_τ, springs = err_springs)
end

const history_rtol = 1.0e-9

# Checks the package does not yet pass, with the audit IDs in
# `docs/derivations/tensor_reduction_audit.md` that explain them. `:solve` is τII
# from the reference history, `:recovery` the spring tensors from the reference
# history, `:loop` both quantities along the package's own time loop.
const known_failures = let
    all_loads = keys(ElasticHistoryCases.load_cases)
    tensor_loads = (:reversal, :rotation, :three_D)
    linear_blocks = (:KV, :mixed, :Burgers, :two_blocks, :spring_and_subbranch, :two_springs_subbranch, :PL_outer)
    nonlinear_blocks = (:PL_direct, :PL_subbranch)
    solve = Set{Tuple{Symbol, Symbol}}()
    recovery = Set{Tuple{Symbol, Symbol}}()
    # S4: power-law viscosities in a block evaluated at the shifted block strain
    # rate. S25: under scalar reversal the block strain rate turns negative,
    # where the scalar power law throws a `DomainError`.
    for t in nonlinear_blocks, l in all_loads
        push!(solve, (t, l))
    end
    # S6/S7: block spring stresses recovered from invariants, coaxial with ε + H.
    for t in linear_blocks, l in tensor_loads
        push!(recovery, (t, l))
    end
    # S5: block power-law viscosities evaluated at ε = 0 during recovery.
    for t in nonlinear_blocks, l in all_loads
        push!(recovery, (t, l))
    end
    (; solve, recovery, loop = union(solve, recovery))
end

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
    check(broken, r) = broken ? (@test_broken r < history_rtol) : (@test r < history_rtol)
    for (tn, c) in pairs(merge(C.linear_topologies, C.nonlinear_topologies)), (ln, load) in pairs(C.load_cases)
        @testset "$tn / $ln" begin
            fed = compare_histories(c, load; feed = :reference)
            own = compare_histories(c, load; feed = :package)
            check((tn, ln) in known_failures.solve, fed.τII)
            check((tn, ln) in known_failures.recovery, fed.springs)
            check((tn, ln) in known_failures.loop, max(own.τII, own.springs))
        end
    end
end
