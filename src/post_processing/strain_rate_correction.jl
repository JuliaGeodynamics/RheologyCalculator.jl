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
    return _history_tensor(c.leafs, c.branches, ε, others.τ0, others, εb)
end

# Split into two independent contributions and add them:
#   1. Direct elastic leafs of the outer SeriesModel (simple Maxwell case).
#   2. Elastic elements inside ParallelModel branches (KV / generalized Maxwell).
# The τ0 tuple is ordered: leafs first (indexed 1 … n_el_leafs), then branches
# (indexed n_el_leafs+1 … end), matching global_eltype_numbering. `εb` holds, per
# branch, the strain rate at which its viscosities are evaluated.
@inline function _history_tensor(leafs::NTuple{N, Any}, branches::NTuple{Nb, Any}, ε, τ0::NTuple{Nτ}, others, εb) where {N, Nb, Nτ}
    n_el_leafs = count_elastic(leafs)
    cor_leafs = if iszero(n_el_leafs)
        ε .* 0
    else
        effective_strain_rate_correction(leafs, (), ε, τ0, others)
    end
    cor_branch = _kv_corrections(branches, ε, τ0, others, n_el_leafs, εb)
    return cor_leafs .+ cor_branch
end

"""
    _block_strain_rates(eqs, x, branches)

The block strain-rate unknowns of `x`, one per branch of the outer
`SeriesModel`. They are the children of its global `compute_strain_rate`
equation that follow its local children, in branch order.
"""
@generated function _block_strain_rates(eqs::NTuple{N, Any}, x, branches::NTuple{Nb, Any}) where {N, Nb}
    k = findfirst(E -> E.parameters[1] === true && E.parameters[3] === typeof(compute_strain_rate), collect(eqs.parameters))
    k === nothing && return :(error("composite has no global deviatoric equation"))
    entries = [:(x[child[end - $Nb + $i]]) for i in 1:Nb]
    return quote
        @inline
        child = eqs[$k].child
        return ($(entries...),)
    end
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

# For an elastic leaf: correction = τ0 / (2η).  This is the Maxwell backstress
# term: the strain the spring would have produced at strain rate τ0/(2G·dt)
# times dt, which must be subtracted from the total ε before solving.
# η = G·dt for an Elasticity element, so τ0/(2η) = τ0/(2G·dt).
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
# Background (see docs/derivations/tensor_reduction.typ for the full derivation):
#
# For a SeriesModel whose branches contain ParallelModel elements, each
# parallel block contributes an effective strain-rate correction beyond the
# simple Maxwell leaf correction already handled by the existing code.
#
# For one ParallelModel branch the corrected strain rate is:
#
#   ε_eff = ε + Σ_i η_star_i * τ0_i / (2 * η_KV)           (*)
#
# where the sum runs over every elastic source (leaf or sub-branch) inside
# the parallel block, and:
#
#   η_KV  = Σ η_eff_i   (sum of effective viscosities, i.e. arithmetic mean)
#   η_star = 1                                  for a direct elastic leaf
#   η_star = η_eff_M / (G * dt) = η_v / (η_v + G*dt)  for a Maxwell sub-branch
#
# η_star < 1 for Maxwell branches: a softer spring (small G) contributes more
# of its backstress; a very stiff spring makes η_star → 0 because the elastic
# strain is negligible and the element behaves like a pure dashpot.
#
# All index arithmetic and η_star expressions are resolved at *compile time*
# by the generated functions below, so the runtime code is a flat sequence of
# arithmetic with no dynamic dispatch or branching.
# -----------------------------------------------------------------------

"""
    count_elastic(r::NTuple{N, AbstractRheology})

Return, at compile time, the number of `AbstractElasticity` elements in the
leaf tuple `r`. Used to compute τ0 index offsets before processing the
`ParallelModel` branches of a `SeriesModel`.
"""
@generated function count_elastic(r::NTuple{N, AbstractRheology}) where {N}
    # Count by inspecting the concrete element types at specialisation time.
    n = count(T -> T <: AbstractElasticity, r.parameters)
    return :($n)
end
count_elastic(::Tuple{}) = 0

"""
    _n_elastic_in_parallel(::Type{ParallelModel{L, B}})

Return, at the *type level*, the total number of elastic elements inside a
`ParallelModel`: elastic direct leafs (type `L`) plus elastic leafs of every
`SeriesModel` sub-branch in `B`.

Called only from `_kv_corrections` at specialisation time; never called at
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

Throw at specialisation time when a branch that carries elastic history has such
a nested composite, whether or not the nested composite is itself elastic.
Branches without elastic elements take no correction and are not restricted.
"""
function _assert_kv_nesting_supported(::Type{ParallelModel{L, B}}) where {L, B}
    _type_has_elastic(ParallelModel{L, B}) || return nothing
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

Per-branch starting offsets into the `τ0` tuple, resolved at specialisation time.

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
    _kv_corrections(branches, ε, τ0, others, offset, εb)

Accumulate the generalized Maxwell / KV effective strain-rate corrections from
all `ParallelModel` branches of a `SeriesModel`.

`offset` is the number of elastic elements already consumed by the series
leafs, so `τ0[offset + k]` is the backstress for the k-th elastic element
inside the branches. Branch `i` evaluates its viscosities at strain rate `εb[i]`.

The τ0 index for each branch is pre-computed at specialisation time (via
`_n_elastic_in_parallel`) and baked in as a literal integer, yielding
allocation-free, branch-free runtime code.
"""
@generated function _kv_corrections(
        branches::NTuple{Nb, Any}, ε, τ0, others, offset, εb
    ) where {Nb}
    offsets = _branch_tau0_offsets(branches)

    # Emit one _kv_branch_correction call per branch, each with its τ0 start
    # index baked in as a literal integer — no runtime bookkeeping needed.
    stmts = Any[:(cor = ε .* 0)]
    for i in 1:Nb
        push!(stmts, :(cor = cor .+ _kv_branch_correction(branches[$i], ε, τ0, others, offset + $(offsets[i]), εb[$i])))
    end
    push!(stmts, :(cor))

    return quote
        @inline
        $(stmts...)
    end
end

"""
    _kv_branch_correction(branch::ParallelModel, ε, τ0, others, el_idx_start, εb)

Compute the generalized Maxwell / Kelvin-Voigt effective strain-rate correction
for a single `ParallelModel` branch.  Returns zero immediately when the branch
contains no elastic elements.

The correction follows equation (*) in `docs/derivations/tensor_reduction.typ`:

    Σ_i η_star_i * τ0_i / (2 * η_KV)

`el_idx_start` is the 1-based index of the elastic element immediately before
the first elastic element owned by this branch (i.e. `τ0[el_idx_start + 1]`
is this branch's first backstress entry). The branch viscosities are evaluated
at the scalar strain rate `εb`.
"""
@inline function _kv_branch_correction(branch::ParallelModel, ε, τ0, others, el_idx_start, εb)
    # Short-circuit: no elastic elements anywhere in this parallel block.
    iselastic(branch) == Val(false) && return ε .* 0
    args = merge((; ε = εb), others)
    # η_KV: arithmetic sum of effective viscosities of all sub-elements
    # (viscous leafs + Maxwell sub-branches).  This is the denominator of (*).
    η_KV = _checked_η_KV(branch.leafs, branch.branches, args)
    # ws = Σ_i η_star_i * τ0_i — the weighted backstress numerator of (*).
    # Same shape as ε (tensor), because τ0 entries are also stored as tensors.
    ws = _weighted_backstress(branch.leafs, branch.branches, ε, τ0, args, el_idx_start)
    # Final correction tensor: ws / (2 * η_KV) broadcasted element-wise.
    return @. ws / (2 * η_KV)
end

"""
    _η_KV(leafs, subs, args)

Compute the Kelvin-Voigt effective viscosity of a parallel block:

    η_KV = Σ_j η_j  (viscous leafs)  +  Σ_k η_eff_M_k  (Maxwell sub-branches)

where `η_eff_M = 1 / (1/η_v + 1/(G*dt))` is the harmonic-mean effective
viscosity of each Maxwell `SeriesModel` sub-branch.
"""
@inline function _η_KV(leafs::NTuple{N, AbstractRheology}, subs::Tuple, args) where {N}
    # Arithmetic sum over the direct leafs of the ParallelModel: a viscous leaf
    # contributes η, an elastic one G*dt. The `false` seed is an additive zero
    # that takes the viscosities' type, so Float32 and dual numbers are preserved.
    η = foldtuple(+, false, l -> compute_viscosity_series(l, args), leafs)
    # Each Maxwell SeriesModel sub-branch contributes its harmonic-mean
    # effective viscosity η_eff_M = (η_v * G*dt)/(η_v + G*dt).
    return foldtuple(+, η, sub -> _η_eff_maxwell(sub.leafs, args), subs)
end

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
otherwise produce a silent `Inf` or `NaN` correction. `Inf` itself is a
legitimate value: a rigid leaf (a plastic element below yield) in parallel with a
spring makes the whole branch rigid, and its correction is correctly zero.
"""
@inline function _checked_η_KV(leafs, subs, args)
    η_KV = _η_KV(leafs, subs, args)
    η_KV > 0 || _throw_nonpositive_η_KV(η_KV)
    return η_KV
end

"""
    _η_eff_maxwell(leafs::NTuple{N, AbstractRheology}, args)

Effective viscosity of a Maxwell `SeriesModel` branch: the harmonic mean of
the viscosities of its constituent leaf elements.

    η_eff_M = 1 / Σ_i (1 / η_i)

For a two-element branch `(viscous, elastic)` this reduces to the standard
Maxwell formula `η_v * G * dt / (η_v + G * dt)`.
"""
# The harmonic mean -- the inverse of the sum of inverses -- is the effective
# viscosity of elements in series.
@inline _η_eff_maxwell(leafs::NTuple{N, AbstractRheology}, args) where {N} =
    inv(foldtuple(+, false, l -> inv(compute_viscosity_series(l, args)), leafs))

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
    _accumulate_weighted_backstress(identity, ε .* 0, leafs, subs, τ0, args, el_idx_start)

# Σ_i η_star_i * entry(τ0_i), accumulated from `seed`. Every τ0 index and every
# η_star expression is resolved here at specialisation time, so the emitted code
# is a flat sequence of multiply-adds with literal indices.
@generated function _accumulate_weighted_backstress(
        entry::E, seed, leafs::NTuple{N, AbstractRheology}, subs::NTuple{Ns, Any}, τ0, args, el_idx_start
    ) where {E, N, Ns}
    direct, per_sub = _elastic_source_positions(leafs, subs)

    stmts = Any[:(ws = seed)]
    el_count = 0

    # Direct elastic leafs of the ParallelModel: η_star = 1, the pure
    # Kelvin-Voigt case, where the backstress enters undiluted.
    for _ in direct
        el_count += 1
        push!(stmts, :(ws = ws .+ entry(τ0[el_idx_start + $el_count])))
    end

    # Springs of a SeriesModel sub-branch: η_star = η_eff_M / (G*dt), with η_eff_M
    # shared by the sub-branch and G the spring's own modulus.
    for j in 1:Ns
        isempty(per_sub[j]) && continue
        η_eff_M = Symbol(:η_eff_M_, j)
        push!(stmts, :($η_eff_M = _η_eff_maxwell(subs[$j].leafs, args)))
        for pos in per_sub[j]
            el_count += 1
            push!(
                stmts,
                :(ws = ws .+ ($η_eff_M / compute_viscosity(subs[$j].leafs[$pos], args)) .* entry(τ0[el_idx_start + $el_count]))
            )
        end
    end

    push!(stmts, :(ws))
    return quote
        @inline
        $(stmts...)
    end
end
