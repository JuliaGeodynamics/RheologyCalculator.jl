"""
    second_invariant(a)
    second_invariant(xx, yy, xy)
    second_invariant(xx, yy, zz, yz, xz, xy)

Return `a` for a scalar invariant, or the second invariant of a 2D or 3D
symmetric deviatoric tensor stored in Voigt-like component order.
"""
@inline second_invariant(a::Number) = a
@inline second_invariant(xx, yy, xy) = √((xx^2 + yy^2 + (-xx - yy)^2) / 2 + xy^2)
@inline second_invariant(xx, yy, zz, yz, xz, xy) = √(0.5 * (xx^2 + yy^2 + zz^2) + xy^2 + yz^2 + xz^2)
# Convenience wrappers: accept either a bare scalar or a Voigt-ordered NTuple.
@inline second_invariant_value(a::Number) = second_invariant(a)
@inline second_invariant_value(a::NTuple) = second_invariant(a...)

# Squared second invariant of a Voigt tuple; the `/ 2` keeps the element type.
@inline _second_invariant_squared(xx, yy, xy) = (xx^2 + yy^2 + (-xx - yy)^2) / 2 + xy^2
@inline _second_invariant_squared(xx, yy, zz, yz, xz, xy) = (xx^2 + yy^2 + zz^2) / 2 + xy^2 + yz^2 + xz^2

"""
    differentiable_second_invariant(a)

`second_invariant_value` with a derivative that stays finite at a zero tensor.

The invariant is not differentiable where the tensor vanishes, which a load
reversal can pass through. There it returns zero with zero derivative instead of
the `NaN` that `√` produces. A scalar or a one-component tuple is returned as a
signed scalar, as `second_invariant_value` does.
"""
@inline differentiable_second_invariant(a::Number) = a
@inline differentiable_second_invariant(a::NTuple{1}) = a[1]
@inline differentiable_second_invariant(a::NTuple) = safe_sqrt(_second_invariant_squared(a...))

# -----------------------------------------------------------------------
# History tensor H
# -----------------------------------------------------------------------
#
# The elastic state functions are τ0-free, so the backstress histories enter the
# solve only through the global deviatoric equation, which is posed on the
# invariant of the assembled tensor
#
#   E = ε + H,   H = Σ_{springs k of the outer series} τ0_k / (2 G_k Δt)
#                  + Σ_{parallel blocks b} Σ_{springs i of b} η*_i τ0_i / (2 η_KV,b)
#
# (see docs/derivations/tensor_reduction.typ and noncoaxial_reduction.typ).
# H has the shape of ε: a Voigt tuple, or a signed scalar for scalar ε.
# -----------------------------------------------------------------------

"""
    effective_strain_rate_correction(c, x, ε, others)
    effective_strain_rate_correction(c, ε, τ0, others)

Return the history tensor `H` of composite `c`, the strain-rate correction that
carries the elastic backstresses into the solve: the global deviatoric equation
of [`solve`](@ref) is posed on the second invariant of `ε + H`.

```
H = Σ_k τ0_k / (2 G_k Δt) + Σ_b Σ_{i ∈ b} η*_i τ0_i / (2 η_KV,b)
```

The first sum runs over the springs of the outer `SeriesModel`, the second over
the springs of each `ParallelModel` block `b`. In a block, `η_KV,b` is the sum of
the effective viscosities of its elements (a sub-branch counts with the inverse
sum of its own), `η*_i = 1` for a direct spring and `η*_i = η_M / (G_i Δt)` for a
spring of a sub-branch with effective viscosity `η_M`. `H` has the shape of `ε`:
a Voigt tuple, or a signed scalar.

The first form reads `τ0` from `others` and evaluates block viscosities at the
block strain-rate unknowns of the solver vector `x`; it applies to every
composite. The second form takes `τ0` explicitly and needs no solver vector, so
it applies only when every element of every block carrying elastic history has
a viscosity that does not depend on the state (see
[`viscosity_depends_on_state`](@ref)); otherwise it throws.
"""
effective_strain_rate_correction(c::SeriesModel, x::AbstractVector, ε, others) =
    _history_tensor(c, generate_equations(c), x, ε, others)

effective_strain_rate_correction(c::SeriesModel, ε, τ0, others) =
    effective_strain_rate_correction(iselastic(c), c, ε, τ0, others)

function effective_strain_rate_correction(::Val{true}, c::SeriesModel, ε, τ0, others)
    _history_coefficients_constant(c) || _throw_history_needs_x()
    # The block viscosities are state-independent, so any strain rate serves.
    εII = second_invariant_value(ε)
    εb = map(_ -> εII, c.branches)
    return _history_tensor(c.leafs, c.branches, ε, τ0, others, εb)
end

# No elastic element anywhere in the composite: return a scalar zero of ε's element
# type, which broadcasts against ε at the call site.
# ε is one entry per tensor component and τ0 one entry per elastic element, so their
# lengths are independent: a purely viscous composite has τ0 = () alongside a
# 3-component ε. Nothing here reads τ0, so it carries no annotation beyond `Tuple`.
@inline effective_strain_rate_correction(::Val{false}, ::SeriesModel, ::NTuple{N, T}, ::Tuple, ::Any) where {N, T} = zero(T)
@inline effective_strain_rate_correction(::Val{false}, ::SeriesModel, ε::Number, ::Tuple, ::Any) = zero(ε)

@noinline function _throw_history_needs_x()
    throw(
        ArgumentError(
            "a parallel block carrying elastic history has an element whose viscosity depends " *
                "on the state, so its history tensor depends on the solution; use " *
                "`effective_strain_rate_correction(c, x, ε, others)` with the solver vector `x`"
        )
    )
end

# True when no block that carries elastic history holds an element whose
# viscosity depends on the state.
@inline _history_coefficients_constant(c::SeriesModel) =
    !foldtuple(|, false, _block_depends_on_state, c.branches)

@inline _block_depends_on_state(b::ParallelModel) =
    _iselastic(b) && (_any_depends_on_state(b.leafs) || foldtuple(|, false, sub -> _any_depends_on_state(sub.leafs), b.branches))

@inline _any_depends_on_state(leafs::Tuple) = foldtuple(|, false, viscosity_depends_on_state, leafs)

# H at solver vector x. Block viscosities are evaluated at the block strain-rate
# unknowns of x.
@inline _history_tensor(c::SeriesModel, eqs, x, ε, others) = _history_tensor(iselastic(c), c, eqs, x, ε, others)
@inline _history_tensor(::Val{false}, ::SeriesModel, eqs, x, ε, others) = ε .* zero(eltype(x))

@inline function _history_tensor(::Val{true}, c::SeriesModel, eqs, x, ε, others)
    # An absent `τ0` is a zero history, as for the element state functions,
    # whose `τ0` and `P0` keywords default to zero.
    hasfield(typeof(others), :τ0) || return ε .* zero(eltype(x))
    εb = _block_strain_rates(eqs, x, c.branches)
    phys = _block_closure_values(eqs, x, c.branches, global_eltype_numbering(c)[2])
    return _history_tensor(c.leafs, c.branches, ε, others.τ0, others, εb, phys)
end

# Split into two independent contributions and add them:
#   1. Direct elastic leafs of the outer SeriesModel (simple Maxwell case).
#   2. Elastic elements inside ParallelModel branches (KV / generalized Maxwell).
# The τ0 tuple is ordered: leafs first (indexed 1 … n_el_leafs), then branches
# (indexed n_el_leafs+1 … end), matching global_eltype_numbering. `εb` holds, per
# branch, the strain rate at which its state-independent viscosities are
# evaluated, and `phys` the closure unknowns of the branch, or `nothing` (see
# `_block_closure_values`).
@inline function _history_tensor(leafs::NTuple{N, Any}, branches::NTuple{Nb, Any}, ε, τ0::NTuple{Nτ}, others, εb, phys = maptuple(_ -> nothing, branches)) where {N, Nb, Nτ}
    n_el_leafs = count_elastic(leafs)
    cor_leafs = if iszero(n_el_leafs)
        ε .* 0
    else
        effective_strain_rate_correction(leafs, (), ε, τ0, others)
    end
    cor_branch = _kv_corrections(branches, ε, τ0, others, n_el_leafs, εb, phys)
    return cor_leafs .+ cor_branch
end

"""
    _block_strain_rates(eqs, x, branches)

The block strain-rate unknowns of `x`, one per branch of the outer
`SeriesModel`. They are the children of its global `compute_strain_rate`
equation that follow its local children, in branch order.
"""
@inline _block_strain_rates(eqs, x, branches) = maptuple(i -> x[i], block_equation_indices(eqs, branches))

"""
    block_equation_indices(eqs, branches)

Indices of the block equations, one per branch of the outer `SeriesModel`: the
children of its global `compute_strain_rate` equation that follow its local
children, in branch order.
"""
@generated function block_equation_indices(eqs::NTuple{N, Any}, branches::NTuple{Nb, Any}) where {N, Nb}
    k = findfirst(E -> E.parameters[1] === true && E.parameters[3] === typeof(compute_strain_rate), collect(eqs.parameters))
    k === nothing && return :(error("composite has no global deviatoric equation"))
    entries = [:(child[end - $Nb + $i]) for i in 1:Nb]
    return quote
        @inline
        child = eqs[$k].child
        return ($(entries...),)
    end
end

"""
    _block_closure_values(eqs, x, branches, el_nums)

Per branch of the outer `SeriesModel`, `nothing` when the branch has no closure
unknowns, otherwise `(; e, t, el)`: the value of `e_b` and of `t_j` per
sub-branch (zero where absent) and the element numbering `el` of the branch.
The closure unknowns follow the block equation, in the order of
`closure_unknowns`.
"""
@generated function _block_closure_values(eqs::NTuple{N, Any}, x, branches::NTuple{Nb, Any}, el_nums) where {N, Nb}
    entries = [
        :(_closure_values(closure_unknowns(branches[$i], compute_stress), branches[$i], iblock[$i], x, el_nums[$i]))
            for i in 1:Nb
    ]
    return quote
        @inline
        iblock = block_equation_indices(eqs, branches)
        return ($(entries...),)
    end
end

@inline _closure_values(::Tuple{}, branch, iblock, x, el) = nothing
@inline function _closure_values(markers::NTuple{N, Any}, branch, iblock, x, el) where {N}
    e_offset, t_offsets = closure_offsets(markers, branch.branches)
    e = _x_at(x, _index_or_zero(_closure_index(e_offset, iblock)))
    t = maptuple(k -> _x_at(x, _index_or_zero(_closure_index(k, iblock))), t_offsets)
    return (; e, t, el)
end

# Scan the leaf tuple for elastic elements and accumulate their corrections.
# The second argument `::Tuple{}` signals "no branches" so only leafs are handled.
# `i` is a running counter incremented only when an elastic leaf is found;
# non-elastic leafs leave i unchanged and contribute nothing (the if-guard
# prevents any τ0 access for them).  Because i is updated by a @nexprs macro
# (resolved at compile time to N literal steps), the loop has no overhead.
@generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, ε, τ0::NTuple{Nτ}, others) where {N, Nτ}
    return quote
        @inline
        i = 0
        ε_elastic_cor = ε .* 0
        Base.@nexprs $N j -> begin
            i = update_correction_index(leafs[j], i)
            if i > 0
                ε_elastic_cor = ε_elastic_cor .+ effective_strain_rate_correction(leafs[j], ε, τ0[i], others, i)
            end
        end

        return ε_elastic_cor
    end
end

# Per-leaf dispatch: route through iselastic so non-elastic elements are a no-op.
@inline effective_strain_rate_correction(c::AbstractRheology, ε, τ0, others, I) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others, I)
@inline effective_strain_rate_correction(::Val{false}, c::AbstractRheology, ε, τ0, others, I) = 0

# A spring of the outer SeriesModel adds τ0 / (2η) = τ0 / (2 G dt) to ε: its
# backward-Euler law ε_s = (τ - τ0) / (2 G dt) moves the backstress term to the
# applied side, so the series sees the strain rate ε + τ0 / (2 G dt).
@inline function effective_strain_rate_correction(::Val{true}, c::AbstractRheology, ε, τ0, others, I)
    η = compute_viscosity(c, merge((; ε), others))
    correction = @. τ0 / (2 * η)
    return correction
end

# Increment the elastic-element counter only when the leaf IS elastic.
# This is used instead of a runtime branch so the @nexprs unrolled loop
# can use a literal τ0 index (`τ0[i]`) without heap allocation.
@inline update_correction_index(c::AbstractRheology, I) = update_correction_index(iselastic(c), I)
@inline update_correction_index(::Val{false}, I) = I       # non-elastic: skip
@inline update_correction_index(::Val{true}, I) = I + 1  # elastic: advance τ0 cursor

# Public wrapper: returns a Val so callers can dispatch on the result without
# paying for a runtime branch (the Val is always resolved at compile time when
# the concrete type of `r` is known, which it always is in @generated contexts).
"""
    iselastic(r)

Return `Val(true)` when `r` is an elastic rheology or a composite containing an
elastic rheology, otherwise `Val(false)`.
"""
@inline iselastic(r::AbstractCompositeModel) = Val(_iselastic(r))
@inline iselastic(::AbstractElasticity) = Val(true)
@inline iselastic(::AbstractRheology) = Val(false)

# Recursive Bool-valued predicate used by the Val-returning wrappers above.
# Checks leafs first, then branches (short-circuits on first true).
@inline _iselastic(r::AbstractCompositeModel) = _iselastic(r.leafs) || _iselastic(r.branches)

# Scan a tuple of leafs or of nested composites (a `branches` field) and return
# true as soon as any entry is or contains an elastic element. @nexprs unrolls
# at compile time; the early return skips the rest.
@inline _iselastic(r::NTuple{N, Union{AbstractRheology, AbstractCompositeModel}}) where {N} =
    _iselastic(first(r)) || _iselastic(Base.tail(r))

@inline _iselastic(::Tuple{}) = false
@inline _iselastic(::AbstractElasticity) = true
@inline _iselastic(::AbstractRheology) = false

# -----------------------------------------------------------------------
# Generalized Maxwell / Kelvin-Voigt strain-rate correction for branches
# -----------------------------------------------------------------------
#
# Derivation: docs/derivations/tensor_reduction.typ (generalized parallel block)
# and noncoaxial_reduction.typ.
#
# Each ParallelModel block b of the outer SeriesModel adds to H the term
#
#   Σ_{springs i of b} η_star_i τ0_i / (2 η_KV)                    (*)
#
# with
#
#   η_KV   = Σ_{direct elements} η + Σ_{sub-branches j} η_M,j   (a sum)
#   η_M,j  = 1 / Σ_{elements of j} (1 / η)                       (an inverse sum)
#   η_star = 1                        for a direct spring
#   η_star = η_M,j / (G_i dt)         for spring i of sub-branch j, with its own G_i
#
# where every element counts with its effective viscosity (a spring with
# G dt). η_star ≤ 1 inside a sub-branch: a stiff spring takes a small share of
# the strain rate, so its backstress weighs little. Elements whose viscosity
# depends on the state are evaluated at the closure unknowns of the block
# (src/equation_system/closure_unknowns.jl).
#
# All index arithmetic and η_star expressions are resolved at compile time by
# the generated functions below, so the runtime code is a flat sequence of
# arithmetic with no dynamic dispatch or branching.
# -----------------------------------------------------------------------

"""
    count_elastic(r::NTuple{N, AbstractRheology})

Return, at compile time, the number of `AbstractElasticity` elements in the
leaf tuple `r`. Used to compute τ0 index offsets before processing the
`ParallelModel` branches of a `SeriesModel`.
"""
@generated function count_elastic(r::NTuple{N, AbstractRheology}) where {N}
    # Count by inspecting the concrete element types at specialization time.
    n = count(T -> T <: AbstractElasticity, r.parameters)
    return :($n)
end
count_elastic(::Tuple{}) = 0

"""
    _n_elastic_in_parallel(::Type{ParallelModel{L, B}})

Return, at the *type level*, the total number of elastic elements inside a
`ParallelModel`: elastic direct leafs (type `L`) plus elastic leafs of every
`SeriesModel` sub-branch in `B`.

Called only from `_kv_corrections` at specialization time; never called at
runtime.
"""
function _n_elastic_in_parallel(::Type{ParallelModel{L, B}}) where {L, B}
    direct, per_sub = _elastic_source_positions(L, B)
    return length(direct) + sum(length, per_sub; init = 0)
end

# Recursively check, at the type level, whether a nested rheology/composite
# type contains any `AbstractElasticity` element anywhere in its subtree.
_type_has_elastic(::Type{<:AbstractElasticity}) = true
_type_has_elastic(::Type{<:AbstractRheology}) = false
function _type_has_elastic(::Type{<:Union{SeriesModel{L, B}, ParallelModel{L, B}}}) where {L, B}
    return any(_type_has_elastic, L.parameters) || any(_type_has_elastic, B.parameters)
end

"""
    _assert_kv_nesting_supported(::Type{ParallelModel{L, B}})

The `_η_KV`/`_η_eff_maxwell`/`_weighted_backstress`/`_n_elastic_in_parallel`
formulas (see `docs/derivations/tensor_reduction.typ`) are derived for a
`ParallelModel` branch whose `SeriesModel` sub-branches contain plain rheology
leafs only, i.e. at most one level of Series/Parallel alternation. `_η_eff_maxwell`
reads only a sub-branch's leafs, so any composite nested inside a sub-branch
would be left out of η_KV and of the backstress weights.

Throw at specialization time when a branch that carries elastic history has such
a nested composite, whether or not the nested composite is itself elastic, or
holds a plastic element: the plastic local equations see the history-free
block strain rate and stress, so the yield condition would ignore the
backstresses. Branches without elastic elements take no correction and are not
restricted.
"""
function _assert_kv_nesting_supported(::Type{ParallelModel{L, B}}) where {L, B}
    _type_has_elastic(ParallelModel{L, B}) || return nothing
    plastic = findfirst(T -> T <: AbstractPlasticity, [L.parameters..., (T for S in B.parameters for T in S.parameters[1].parameters)...])
    plastic === nothing || error(
        "elastic strain-rate correction: a ParallelModel that carries elastic history holds a " *
            "plastic element. Plastic elements are not supported in a ParallelModel with a spring; " *
            "place the plastic element in the outer SeriesModel or in a block without springs."
    )
    for S in B.parameters
        sub_branches = S.parameters[2]  # SeriesModel sub-branch's own `branches` field type
        isempty(sub_branches.parameters) || error(
            "elastic strain-rate correction: a ParallelModel that carries elastic history has a " *
                "SeriesModel sub-branch with a nested composite ($(first(sub_branches.parameters))). " *
                "Sub-branches of such a ParallelModel may contain rheology elements only " *
                "(see docs/derivations/tensor_reduction.typ)."
        )
    end
    return nothing
end

"""
    _branch_tau0_offsets(branches::Type)

Per-branch starting offsets into the `τ0` tuple, resolved at specialization time.

`τ0` is ordered to match `global_eltype_numbering`: the outer `SeriesModel`'s own
elastic leafs first, then its branches left to right, each branch owning
`_n_elastic_in_parallel` consecutive entries. So branch `i` reads its k-th
backstress from `τ0[offset + offsets[i] + k]`, where `offset` is what the series
leafs consumed.

Also checks that no branch nests an elastic element deeper than the correction
formulas are derived for.
"""
function _branch_tau0_offsets(branches::Type)
    foreach(_assert_kv_nesting_supported, branches.parameters)
    counts = [_n_elastic_in_parallel(T) for T in branches.parameters]
    return cumsum([0; counts[1:(end - 1)]])
end

"""
    _kv_corrections(branches, ε, τ0, others, offset, εb[, phys])

Accumulate the generalized Maxwell / KV effective strain-rate corrections from
all `ParallelModel` branches of a `SeriesModel`.

`offset` is the number of elastic elements already consumed by the series
leafs, so `τ0[offset + k]` is the backstress for the k-th elastic element
inside the branches. Branch `i` evaluates its state-independent viscosities at
strain rate `εb[i]`, and the others at its closure unknowns `phys[i]` (see
`_block_closure_values`); without `phys` every viscosity is evaluated at `εb`.

The τ0 index for each branch is pre-computed at specialization time (via
`_n_elastic_in_parallel`) and baked in as a literal integer, yielding
allocation-free, branch-free runtime code.
"""
@inline _kv_corrections(branches::Tuple, ε, τ0, others, offset, εb) =
    _kv_corrections(branches, ε, τ0, others, offset, εb, maptuple(_ -> nothing, branches))

@generated function _kv_corrections(
        branches::NTuple{Nb, Any}, ε, τ0, others, offset, εb, phys
    ) where {Nb}
    offsets = _branch_tau0_offsets(branches)

    # Emit one _kv_branch_correction call per branch, each with its τ0 start
    # index baked in as a literal integer — no runtime bookkeeping needed.
    stmts = Any[:(cor = ε .* 0)]
    for i in 1:Nb
        push!(stmts, :(cor = cor .+ _kv_branch_correction(branches[$i], ε, τ0, others, offset + $(offsets[i]), εb[$i], phys[$i])))
    end
    push!(stmts, :(cor))

    return quote
        @inline
        $(stmts...)
    end
end

"""
    _kv_branch_correction(branch::ParallelModel, ε, τ0, others, el_idx_start, εb[, phys])

Compute the generalized Maxwell / Kelvin-Voigt effective strain-rate correction
for a single `ParallelModel` branch.  Returns zero immediately when the branch
contains no elastic elements.

The correction follows equation (*) in `docs/derivations/tensor_reduction.typ`:

    Σ_i η_star_i * τ0_i / (2 * η_KV)

`el_idx_start` is the number of `τ0` entries that precede this branch, so
`τ0[el_idx_start + 1]` is its first backstress. The branch viscosities are evaluated
as in `_block_state`, at the scalar strain rate `εb` and the closure unknowns
`phys`.
"""
@inline function _kv_branch_correction(branch::ParallelModel, ε, τ0, others, el_idx_start, εb, phys = nothing)
    # Short-circuit: no elastic elements anywhere in this parallel block.
    iselastic(branch) == Val(false) && return ε .* 0
    args = merge((; ε = εb), others)
    # η_KV: arithmetic sum of effective viscosities of all sub-elements
    # (viscous leafs + Maxwell sub-branches), the denominator of (*).
    # ws = Σ_i η_star_i * τ0_i, the weighted backstress numerator of (*), with
    # the shape of ε (τ0 entries are stored as tensors).
    η_KV, ws = _block_state(branch.leafs, branch.branches, τ0, args, el_idx_start, phys, ε .* 0)
    return @. ws / (2 * η_KV)
end

"""
    _η_KV(leafs, subs, args[, phys])

Compute the Kelvin-Voigt effective viscosity of a parallel block:

    η_KV = Σ_j η_j  (direct leafs)  +  Σ_k η_eff_M_k  (Maxwell sub-branches)

where a direct dashpot counts with `η`, a direct spring with `G dt`, and
`η_eff_M = 1 / Σ_i (1/η_i)` is the inverse sum over the elements of each
`SeriesModel` sub-branch. Elements are evaluated at `args` and `phys` as
described above `_secant_viscosity`.
"""
# The `false` seed is an additive zero that takes the viscosities' type, so
# Float32 and dual numbers are preserved.
@inline _η_KV(leafs::NTuple{N, AbstractRheology}, subs::Tuple, args, phys = nothing) where {N} =
    _sum_viscosities(_direct_viscosities(leafs, args, phys), _subbranch_viscosities(subs, args, phys))

@inline _sum_viscosities(direct, sub) = foldtuple(+, foldtuple(+, false, identity, direct), identity, sub)

@noinline function _throw_nonpositive_η_KV(η_KV)
    throw(
        ArgumentError(
            "effective Kelvin-Voigt viscosity of this parallel branch is $η_KV; every " *
                "element in a branch carrying elastic history must define compute_viscosity, " *
                "and `others` must supply a nonzero `dt`"
        )
    )
end

"""
    _checked_η_KV(leafs, subs, args)

`_η_KV`, guarded against an aggregate that cannot describe a real branch.

Every caller reaches this only after establishing that the branch carries
elastic history, so a zero, negative, or `NaN` aggregate is a defect — a leaf
with no `compute_viscosity` method, or an `others` that omits `dt` — and would
otherwise produce a silent `Inf` or `NaN` correction. `Inf` passes: an element of
infinite viscosity makes the whole branch rigid, and its correction is zero.
"""
@inline _checked_η_KV(leafs, subs, args, phys = nothing) = _check_η_KV(_η_KV(leafs, subs, args, phys))

@inline function _check_η_KV(η_KV)
    η_KV > 0 || _throw_nonpositive_η_KV(η_KV)
    return η_KV
end

"""
    _η_eff_maxwell(leafs::NTuple{N, AbstractRheology}, args)

Effective viscosity of a Maxwell `SeriesModel` branch: the inverse of the sum
of the inverse viscosities of its leaf elements.

    η_eff_M = 1 / Σ_i (1 / η_i)

For a two-element branch `(viscous, elastic)` this reduces to the standard
Maxwell formula `η_v * G * dt / (η_v + G * dt)`.
"""
# Elements in series share the stress and add their strain rates, so their
# inverse viscosities add.
@inline _η_eff_maxwell(leafs::NTuple{N, AbstractRheology}, args) where {N} =
    _η_eff_maxwell(leafs, args, nothing, nothing)

@inline _η_eff_maxwell(leafs::NTuple{N, AbstractRheology}, args, ::Nothing, ::Nothing) where {N} =
    inv(foldtuple(+, false, l -> inv(compute_viscosity_series(l, args)), leafs))

@inline _η_eff_maxwell(leafs::NTuple{N, AbstractRheology}, args, t, el) where {N} =
    inv(foldtuple(+, false, inv, maptuple((l, n) -> _element_viscosity(compute_strain_rate, l, args, t, n), leafs, el)))

# -----------------------------------------------------------------------
# Element viscosities of a block carrying elastic history
# -----------------------------------------------------------------------
#
# `phys` is `nothing` when every element is evaluated at `args`, whose `ε` is
# the shifted block strain rate; this is exact for elements whose viscosity
# does not depend on the state. Otherwise `phys = (; e, t, el)` holds the
# closure unknowns of the block (see src/equation_system/closure_unknowns.jl)
# and its element numbering, and every element for which `_uses_secant` holds
# is evaluated at its physical argument: `e` for a direct element, `t[j]` for an
# element of sub-branch `j`. The graph equations evaluate these elements the
# same way (`_evaluate_leaf`), so both share one set of coefficients.
# -----------------------------------------------------------------------

"""
    _secant_viscosity(fn, r, args, s, n)

Secant viscosity of element `r` at physical argument `s`, from its own state
function evaluated at `a = |s|`: `compute_stress(r; ε = a) / (2a)` for a
strain-rate argument (a direct element of a block), `a / (2
compute_strain_rate(r; τ = a))` for a stress argument (an element of a
sub-branch). `s` is signed for a scalar strain rate. `n` is the element number
used to extract its history fields from `args`.
"""
@inline function _secant_viscosity(::typeof(compute_stress), r, args, s, n)
    e = abs(s)
    local_args = merge(extract_local_kwargs(args, history_kwargs(r), n), (; ε = e))
    return compute_stress(r, local_args) * safe_inv(2 * e)
end

@inline function _secant_viscosity(::typeof(compute_strain_rate), r, args, s, n)
    t = abs(s)
    local_args = merge(extract_local_kwargs(args, history_kwargs(r), n), (; τ = t))
    return t * safe_inv(2 * compute_strain_rate(r, local_args))
end

@inline _element_viscosity(fn::F, r, args, s, n) where {F} = _element_viscosity(Val(_uses_secant(r)), fn, r, args, s, n)
@inline _element_viscosity(::Val{false}, ::F, r, args, s, n) where {F} = compute_viscosity_series(r, args)
@inline _element_viscosity(::Val{true}, fn::F, r, args, s, n) where {F} = _secant_viscosity(fn, r, args, s, n)

@inline _direct_viscosities(leafs, args, ::Nothing) = maptuple(l -> compute_viscosity_series(l, args), leafs)
@inline _direct_viscosities(leafs, args, phys) =
    maptuple((l, n) -> _element_viscosity(compute_stress, l, args, phys.e, n), leafs, phys.el[1])

@inline _subbranch_physical(::Nothing, j) = (nothing, nothing)
@inline _subbranch_physical(phys, j) = (phys.t[j], phys.el[2][j][1])

# Effective viscosity of each sub-branch, in sub-branch order.
@generated function _subbranch_viscosities(subs::NTuple{Ns, Any}, args, phys) where {Ns}
    entries = [:(_η_eff_maxwell(subs[$j].leafs, args, _subbranch_physical(phys, $j)...)) for j in 1:Ns]
    return quote
        @inline
        $(Expr(:tuple, entries...))
    end
end

# Compile-time positions of the elastic sources inside a `ParallelModel` branch:
# the elastic direct leafs, and per sub-branch the elastic leafs within it.
# `leafs` and `subs` are the *types* seen by a @generated function.
#
# The order is the τ0 order for the branch: direct leafs first, then sub-branches
# left to right, so the k-th source found here owns τ0[el_idx_start + k].
function _elastic_source_positions(leafs::Type, subs::Type)
    direct = findall(Ti -> Ti <: AbstractElasticity, collect(leafs.parameters))
    per_sub = [
        findall(Ti -> Ti <: AbstractElasticity, collect(subs.parameters[j].parameters[1].parameters))
            for j in 1:length(subs.parameters)
    ]
    return direct, per_sub
end

"""
    _weighted_backstress(leafs, subs, ε, τ0, args, el_idx_start)

Compute the weighted backstress numerator `Σ_i η_star_i * τ0_i` for a
`ParallelModel` branch, with the shape of `ε`.

Two classes of elastic sources contribute:
- **Direct elastic leafs** of the `ParallelModel`: `η_star = 1`.  These
  correspond to the simple Kelvin-Voigt case where the elastic element is a
  direct parallel element (backstress enters undiluted).
- **Springs of a `SeriesModel` sub-branch**: `η_star = η_eff_M / η_el`, where
  `η_eff_M` is the effective viscosity of the whole sub-branch and
  `η_el = G * dt` uses that spring's own `G`. This is the generalized Maxwell
  weighting: a softer spring (small G) makes `η_star → 1`; a very stiff spring
  makes `η_star → 0`.

All τ0 index literals and η_star computations are resolved at *compile time*
(via `@generated`), so the emitted code is a flat sequence of multiply-adds.
`el_idx_start` carries the τ0 offset inherited from the outer `SeriesModel`
leaf count.
"""
@inline _weighted_backstress(leafs, subs, ε, τ0, args, el_idx_start) =
    _block_state(leafs, subs, τ0, args, el_idx_start, nothing, ε .* 0)[2]

"""
    _block_state(leafs, subs, τ0, args, el_idx_start, phys, seed)

The coefficients of a `ParallelModel` block with direct elements `leafs` and
sub-branches `subs`, as `(η_KV, W, ηM, S)`:

- `η_KV`: the effective viscosity of the block, checked to be positive;
- `W = Σ_i η_star_i τ0_i`: the weighted backstress of all its springs,
  accumulated from `seed`, which sets its shape;
- `ηM[j]`: the effective viscosity of sub-branch `j`;
- `S[j]`: the part of `W` owned by the springs of sub-branch `j`.

Element viscosities are evaluated at `args` and `phys` as described above
`_secant_viscosity`. Every τ0 index and every η_star expression is resolved at
specialization time, so the emitted code is a flat sequence of multiply-adds
with literal indices.
"""
@generated function _block_state(
        leafs::NTuple{N, AbstractRheology}, subs::NTuple{Ns, Any}, τ0, args, el_idx_start, phys, seed
    ) where {N, Ns}
    direct, per_sub = _elastic_source_positions(leafs, subs)

    stmts = Any[
        :(ηM = _subbranch_viscosities(subs, args, phys)),
        :(η_KV = _check_η_KV(_sum_viscosities(_direct_viscosities(leafs, args, phys), ηM))),
        :(W = seed),
    ]
    el_count = 0

    # Direct elastic leafs of the ParallelModel: η_star = 1, the pure
    # Kelvin-Voigt case, where the backstress enters undiluted.
    for _ in direct
        el_count += 1
        push!(stmts, :(W = W .+ τ0[el_idx_start + $el_count]))
    end

    # Springs of a SeriesModel sub-branch: η_star = η_eff_M / (G*dt), with η_eff_M
    # shared by the sub-branch and G the spring's own modulus.
    S = [Symbol(:S_, j) for j in 1:Ns]
    for j in 1:Ns
        push!(stmts, :($(S[j]) = seed))
        for pos in per_sub[j]
            el_count += 1
            push!(
                stmts,
                :($(S[j]) = $(S[j]) .+ (ηM[$j] / compute_viscosity(subs[$j].leafs[$pos], args)) .* τ0[el_idx_start + $el_count])
            )
        end
        push!(stmts, :(W = W .+ $(S[j])))
    end

    push!(stmts, :((η_KV, W, ηM, $(Expr(:tuple, S...)))))
    return quote
        @inline
        $(stmts...)
    end
end
