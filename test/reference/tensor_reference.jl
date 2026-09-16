# Monolithic full-tensor reference solver for elastic history tests.
#
# One backward-Euler step of a composite `SeriesModel` is solved with every
# unknown held as a full deviatoric tensor, so no invariant reduction or
# strain-rate correction is involved:
#
#   unknowns   τ (outer series stress), ε_b (strain rate of each parallel
#              block b), τ_j (stress of each series sub-branch j of a block)
#   residuals  ε  - Σ_{outer leafs} ε_l(τ)       - Σ_b ε_b   = 0
#              τ  - Σ_{direct leafs of b} τ_l(ε_b) - Σ_j τ_j = 0   (per block)
#              ε_b - Σ_{leafs of j} ε_l(τ_j)                = 0   (per sub-branch)
#
# Element laws are isotropic and built from the package's scalar state
# functions. Springs are linear, so their tensor laws carry the history
# explicitly: ε = (τ - τᵒ)/(2GΔt) in series and τ = τᵒ + 2GΔt ε in parallel.
# Every other element is coaxial with its argument: the scalar law evaluated at
# the second invariant, scaled onto the tensor direction.
#
# Supported topologies: a `SeriesModel` whose branches are `ParallelModel`s,
# whose branches are `SeriesModel`s of leaf elements only. Viscous and elastic
# elements only. Anything else throws.
#
# Tensors are tuples in Voigt order: `(xx, yy, xy)` in 2D, `(xx, yy, zz, yz,
# xz, xy)` in 3D, or a 1-tuple for a signed scalar strain rate.
# The spring history is one tensor per elastic element, in the element order of
# `global_eltype_numbering`: outer leafs, then per block its direct leafs
# followed by the leafs of its sub-branches.
module TensorReference

using ForwardDiff, LinearAlgebra
using RheologyCalculator
using RheologyCalculator: AbstractElasticity, AbstractViscosity,
    SeriesModel, ParallelModel, compute_strain_rate, compute_stress, second_invariant

export ReferenceModel, reference_step, initial_guess, tensor_invariant

"""
    tensor_invariant(t)

Second invariant of a Voigt tuple; the absolute value for a 1-tuple.
"""
tensor_invariant(t::NTuple{1}) = abs(t[1])
tensor_invariant(t::NTuple{3}) = second_invariant(t...)
tensor_invariant(t::NTuple{6}) = second_invariant(t...)

struct ReferenceModel
    outer::Vector{Any}                      # outer series leafs
    blocks::Vector{Vector{Any}}             # direct leafs of each block
    subbranches::Vector{Vector{Vector{Any}}} # leafs of each sub-branch, per block
end

function ReferenceModel(c::SeriesModel)
    blocks = Vector{Any}[]
    subbranches = Vector{Vector{Any}}[]
    for b in c.branches
        b isa ParallelModel || error("reference: outer branches must be ParallelModel, got $(typeof(b))")
        push!(blocks, _checked_leafs(b.leafs))
        subs = Vector{Any}[]
        for j in b.branches
            j isa SeriesModel || error("reference: block branches must be SeriesModel, got $(typeof(j))")
            isempty(j.branches) || error("reference: sub-branches must contain leaf elements only")
            push!(subs, _checked_leafs(j.leafs))
        end
        push!(subbranches, subs)
    end
    return ReferenceModel(_checked_leafs(c.leafs), blocks, subbranches)
end

function _checked_leafs(leafs)
    for r in leafs
        r isa Union{AbstractElasticity, AbstractViscosity} ||
            error("reference: unsupported element $(typeof(r))")
    end
    return collect(Any, leafs)
end

n_springs(leafs) = count(r -> r isa AbstractElasticity, leafs)

function n_springs(m::ReferenceModel)
    n = n_springs(m.outer)
    for (b, subs) in zip(m.blocks, m.subbranches)
        n += n_springs(b) + sum(n_springs, subs; init = 0)
    end
    return n
end

_scale(a, t) = map(ti -> a * ti, t)
_add(a, b) = map(+, a, b)
_sub(a, b) = map(-, a, b)

function _direction_ratio(s, what)
    iszero(s) && error("reference: zero $what invariant; the coaxial law has no direction")
    return s
end

# Strain rate of a series element at stress tensor τ.
series_strain_rate(r::AbstractElasticity, τ, τᵒ, dt) =
    _scale(compute_strain_rate(r; τ = one(eltype(τ)), dt), _sub(τ, τᵒ))
function series_strain_rate(r, τ, τᵒ, dt)
    s = _direction_ratio(tensor_invariant(τ), "stress")
    return _scale(compute_strain_rate(r; τ = s, dt) / s, τ)
end

# Stress of a parallel element at strain-rate tensor ε.
parallel_stress(r::AbstractElasticity, ε, τᵒ, dt) =
    _add(τᵒ, _scale(compute_stress(r; ε = one(eltype(ε)), dt), ε))
function parallel_stress(r, ε, τᵒ, dt)
    e = _direction_ratio(tensor_invariant(ε), "strain-rate")
    return _scale(compute_stress(r; ε = e, dt) / e, ε)
end

# Sum of series strain rates of `leafs` at τ; consumes spring histories from `hist`.
function _series_sum(leafs, τ, hist, k, dt)
    acc = _scale(zero(eltype(τ)), τ)
    for r in leafs
        τᵒ = r isa AbstractElasticity ? hist[k += 1] : nothing
        acc = _add(acc, series_strain_rate(r, τ, τᵒ, dt))
    end
    return acc, k
end

function _parallel_sum(leafs, ε, hist, k, dt)
    acc = _scale(zero(eltype(ε)), ε)
    for r in leafs
        τᵒ = r isa AbstractElasticity ? hist[k += 1] : nothing
        acc = _add(acc, parallel_stress(r, ε, τᵒ, dt))
    end
    return acc, k
end

n_unknown_tensors(m::ReferenceModel) = 1 + length(m.blocks) + sum(length, m.subbranches; init = 0)

_tensor(y, i, nc) = ntuple(c -> y[(i - 1) * nc + c], nc)

# Residual in scaled variables: stresses in units of τc, strain rates of εc.
function residual(m::ReferenceModel, y, ε, hist, dt, τc, εc)
    nc = length(ε)
    T = eltype(y)
    out = Vector{T}(undef, length(y))
    store!(i, t) = for c in 1:nc
        out[(i - 1) * nc + c] = t[c]
    end
    τ = _scale(τc, _tensor(y, 1, nc))
    rtop, k = _series_sum(m.outer, τ, hist, 0, dt)
    rtop = _sub(ε, rtop)
    i = 1
    for (leafs, subs) in zip(m.blocks, m.subbranches)
        εb = _scale(εc, _tensor(y, i += 1, nc))
        ib = i
        rtop = _sub(rtop, εb)
        rb, k = _parallel_sum(leafs, εb, hist, k, dt)
        rb = _sub(τ, rb)
        for jleafs in subs
            τj = _scale(τc, _tensor(y, i += 1, nc))
            rb = _sub(rb, τj)
            rj, k = _series_sum(jleafs, τj, hist, k, dt)
            store!(i, _scale(inv(εc), _sub(εb, rj)))
        end
        store!(ib, _scale(inv(τc), rb))
    end
    store!(1, _scale(inv(εc), rtop))
    return out
end

"""
    reference_step(m, ε, hist, dt; guess, τc, εc) -> (τ, springs, y)

Solve one step for strain-rate tensor `ε` and spring history `hist` (one tensor
per spring). Returns the outer stress tensor, the updated spring tensors in
`global_eltype_numbering` order, and the scaled unknown vector (reusable as the
next `guess`). Throws if Newton does not converge.
"""
function reference_step(m::ReferenceModel, ε::NTuple{nc}, hist, dt; guess, τc, εc) where {nc}
    length(hist) == n_springs(m) || error("reference: expected $(n_springs(m)) spring tensors, got $(length(hist))")
    f(y) = residual(m, y, ε, hist, dt, τc, εc)
    y = collect(float.(guess))
    length(y) == n_unknown_tensors(m) * nc || error("reference: guess has the wrong length")
    converged = false
    for _ in 1:100
        r = f(y)
        J = ForwardDiff.jacobian(f, y)
        Δy = J \ r
        y .-= Δy
        if norm(Δy, Inf) ≤ 1.0e-14 * max(1.0, norm(y, Inf))
            converged = norm(f(y), Inf) ≤ 1.0e-10
            break
        end
    end
    converged || error("reference: Newton did not converge")
    return _recover(m, y, ε, hist, dt, τc, εc)..., y
end

function _recover(m::ReferenceModel, y, ε::NTuple{nc}, hist, dt, τc, εc) where {nc}
    τ = _scale(τc, _tensor(y, 1, nc))
    springs = Any[]
    k = 0
    for r in m.outer
        r isa AbstractElasticity && (k += 1; push!(springs, τ))
    end
    i = 1
    for (leafs, subs) in zip(m.blocks, m.subbranches)
        εb = _scale(εc, _tensor(y, i += 1, nc))
        for r in leafs
            r isa AbstractElasticity && (k += 1; push!(springs, parallel_stress(r, εb, hist[k], dt)))
        end
        for jleafs in subs
            τj = _scale(τc, _tensor(y, i += 1, nc))
            for r in jleafs
                r isa AbstractElasticity && (k += 1; push!(springs, τj))
            end
        end
    end
    return τ, Tuple(springs)
end

"""
    initial_guess(m, ε)

A scaled guess with nonzero invariants everywhere: the outer stress coaxial with
`ε`, each block taking a share of the strain rate and each sub-branch a share of
the stress.
"""
function initial_guess(m::ReferenceModel, ε::NTuple{nc}) where {nc}
    d = _scale(inv(tensor_invariant(ε)), ε)
    nb = length(m.blocks)
    y = Float64[]
    append!(y, d)
    for subs in m.subbranches
        append!(y, _scale(1 / (nb + 1), d))
        for _ in subs
            append!(y, _scale(1 / (length(subs) + 1), d))
        end
    end
    return y
end

end # module
