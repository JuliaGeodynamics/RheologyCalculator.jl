# -----------------------------------------------------------------------
# Closure unknowns of parallel blocks that carry elastic history
# -----------------------------------------------------------------------
#
# Inside a `ParallelModel` block that holds a spring, the graph unknowns are
# shifted values: the block strain rate x̃_b and the sub-branch stresses x̃_j of
# the history-free problem with the same outer stress. With the direction
# n = E / E_II of the assembled strain rate E = ε + H, the physical tensors are
#
#   ε̇_b = x̃_b n - W / (2 η_KV)
#   τ_j = x̃_j n - (η_M,j / η_KV) W + S_j
#
# where W = Σ_{springs i of b} η*_i τ0_i, S_j is the part of W owned by the
# springs of sub-branch j, η_M,j the effective viscosity of sub-branch j and
# η_KV the effective viscosity of the block.
#
# An element whose viscosity depends on the state (see `_uses_secant`) must be
# evaluated at its physical argument, which depends on those same viscosities.
# The loop is closed with physical invariants as extra unknowns:
#
#   e_b = ‖ε̇_b‖   for a block with such an element among its direct leafs
#   t_j = ‖τ_j‖   for a sub-branch with such an element
#
# where ‖·‖ is the second invariant of a tensor and the signed value itself for
# a scalar strain rate, as everywhere else in the solve. Every such element is
# evaluated through its own state function at |e_b| or |t_j|, as a secant
# viscosity (see `_secant_viscosity`). A signed scalar closure unknown lets the
# Newton iteration step across zero after a load reversal; an absolute value
# would have to approach zero, where the secant viscosity of a power law is
# singular. The graph
# equations use the same secant viscosities on the shifted unknowns, so the
# graph, the history tensor H and the closures share one set of coefficients
# and the converged solution satisfies the tensor equations exactly.
#
# The closure unknowns are the first local equations of their block, e_b first,
# then t_j in sub-branch order. Blocks whose elements all have state-independent
# viscosities get none, so their equation layout is unchanged.
# -----------------------------------------------------------------------

"""
    compute_block_strain_rate_invariant

Equation kind of the closure unknown `e_b`: the second invariant of the physical
strain rate of a parallel block that carries elastic history and has a direct
element whose viscosity depends on the state. Its key in [`x_keys`](@ref) is
`:εII_block`.
"""
function compute_block_strain_rate_invariant end

"""
    compute_subbranch_stress_invariant

Equation kind of the closure unknown `t_j`: the second invariant of the physical
stress of a series sub-branch, inside a parallel block that carries elastic
history, that holds an element whose viscosity depends on the state. Its key in
[`x_keys`](@ref) is `:τII_subbranch`.
"""
function compute_subbranch_stress_invariant end

const ClosureFunction = Union{typeof(compute_block_strain_rate_invariant), typeof(compute_subbranch_stress_invariant)}

@inline differentiable_kwargs(::Type{T}, ::typeof(compute_block_strain_rate_invariant)) where {T} = (; εII_block = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_subbranch_stress_invariant)) where {T} = (; τII_subbranch = zero(T))

"""
    _uses_secant(r)

Whether element `r` is evaluated at its physical argument when it sits in a
parallel block that carries elastic history. Springs and plastic elements are
not: springs are linear, and plastic elements are handled by their own local
equations.
"""
@inline _uses_secant(r::AbstractRheology) = viscosity_depends_on_state(r)
@inline _uses_secant(::AbstractElasticity) = false
@inline _uses_secant(::AbstractPlasticity) = false
@inline _any_secant(leafs::Tuple) = foldtuple(|, false, _uses_secant, leafs)

# Markers for the closure unknowns of a block, in equation order.
struct BlockStrainRateClosure end
struct SubbranchStressClosure{J} end

@inline closure_function(::BlockStrainRateClosure) = compute_block_strain_rate_invariant
@inline closure_function(::SubbranchStressClosure) = compute_subbranch_stress_invariant

@inline _closure_subbranch(::BlockStrainRateClosure) = 0
@inline _closure_subbranch(::SubbranchStressClosure{J}) where {J} = J

"""
    closure_unknowns(c, fn)

The closure unknowns that composite `c` adds in the pass of global state
function `fn`, as a tuple of markers. Only a `ParallelModel` that carries
elastic history adds any, and only in the deviatoric pass: one
`BlockStrainRateClosure` when a direct leaf uses a secant viscosity, and one
`SubbranchStressClosure{j}` for each sub-branch `j` with such an element.
"""
@inline closure_unknowns(::AbstractCompositeModel, ::F) where {F} = ()
@inline closure_unknowns(c::ParallelModel, ::typeof(compute_stress)) = _closure_unknowns(iselastic(c), c)
@inline _closure_unknowns(::Val{false}, ::ParallelModel) = ()
@inline _closure_unknowns(::Val{true}, c::ParallelModel) =
    (_block_closure_marker(Val(_any_secant(c.leafs)))..., _subbranch_closure_markers(c.branches)...)

@inline _block_closure_marker(::Val{true}) = (BlockStrainRateClosure(),)
@inline _block_closure_marker(::Val{false}) = ()

@generated function _subbranch_closure_markers(subs::NTuple{Ns, Any}) where {Ns}
    markers = [:(_subbranch_closure_marker(Val(_any_secant(subs[$j].leafs)), Val($j))...) for j in 1:Ns]
    return quote
        @inline
        $(Expr(:tuple, markers...))
    end
end
@inline _subbranch_closure_marker(::Val{true}, ::Val{J}) where {J} = (SubbranchStressClosure{J}(),)
@inline _subbranch_closure_marker(::Val{false}, ::Val) = ()

"""
    closure_offsets(markers, branches)

Positions of the closure unknowns relative to their block equation: `e_offset`
for `e_b` and one entry per branch for `t_j`, each an `Int` or `nothing` when
the unknown is absent. Resolved from the marker types, so whether an unknown
exists is known at compile time.
"""
@generated function closure_offsets(::M, ::NTuple{Nb, Any}) where {M <: Tuple, Nb}
    e = nothing
    t = Any[nothing for _ in 1:Nb]
    for (k, T) in enumerate(M.parameters)
        if T === BlockStrainRateClosure
            e = k
        else
            t[T.parameters[1]] = k
        end
    end
    return :(($e, $(Expr(:tuple, t...))))
end

@inline _closure_index(::Nothing, iblock) = nothing
@inline _closure_index(k::Int, iblock) = iblock + k
@inline _own_closure(::Nothing, iblock, closure) = closure
@inline _own_closure(k::Int, iblock, closure) = iblock + k
@inline _index_or_zero(::Nothing) = 0
@inline _index_or_zero(i::Int) = i

# Entry `i` of `x`, or zero for an absent unknown (`i == 0`).
@inline _x_at(x, i::Int) = iszero(i) ? zero(eltype(x)) : x[i]

"""
    BlockClosure{J, Ns}

Indices into the solver vector that a closure equation of a parallel block
reads: the block strain rate `block`, `e_b` (`0` if absent), the sub-branch
stresses `subs`, and `t_j` per sub-branch (`0` if absent). `τ0_start` is the
`τ0` index just before the block's first spring. `J` is the sub-branch whose
`t_j` the equation solves for, `0` for `e_b`.
"""
struct BlockClosure{J, Ns}
    block::Int
    e::Int
    subs::NTuple{Ns, Int}
    t::NTuple{Ns, Int}
    τ0_start::Int
end

@inline BlockClosure(::Val{J}, block, e, subs::NTuple{Ns, Int}, t::NTuple{Ns, Int}, τ0_start) where {J, Ns} =
    BlockClosure{J, Ns}(block, e, subs, t, τ0_start)

@inline add_closure_equations(::Tuple{}, c, iblock, e_offset, isubs, sub_closures, iself_ref, el_num) = ()

# Left to right, so each closure equation takes the next number from iself_ref.
@inline function add_closure_equations(markers::Tuple, c, iblock, e_offset, isubs, sub_closures, iself_ref, el_num)
    e = _index_or_zero(_closure_index(e_offset, iblock))
    t = maptuple(_index_or_zero, sub_closures)
    τ0_start = _tau0_start(c, el_num)
    return maptuple(markers) do m
        iself_ref[] += 1
        J = _closure_subbranch(m)
        closure = BlockClosure(Val(J), iblock, e, isubs, t, τ0_start)
        CompositeEquation(iblock, (), iself_ref[], closure_function(m), c, J, Val(false), el_num, closure)
    end
end

# The τ0 index just before the first spring of block `c`: springs are numbered
# in `global_eltype_numbering` order, direct leafs first, then sub-branches.
@generated function _tau0_start(::ParallelModel{L, B}, el_num) where {L, B}
    direct, per_sub = _elastic_source_positions(L, B)
    isempty(direct) || return :(el_num[1][$(first(direct))] - 1)
    for j in eachindex(per_sub)
        isempty(per_sub[j]) || return :(el_num[2][$j][1][$(first(per_sub[j]))] - 1)
    end
    return :(error("parallel block has no elastic element"))
end

# A closure equation has no parent or children in the graph, and its rheology
# (the whole block) contributes no state-function term: its residual is added
# by `add_closure_residuals`.
add_child(::SVector{N, T}, ::CompositeEquation{A, B, F}, eq_ind) where {N, T, A, B, F <: ClosureFunction} = zero(T)
@inline subtract_parent(::SVector, ::CompositeEquation{false, T, F}, ::NamedTuple) where {T, F <: ClosureFunction} = false
@inline evaluate_state_function(::ClosureFunction, ::ParallelModel, args, others, el_number) = false

# -----------------------------------------------------------------------
# Secant evaluation in the graph equations
# -----------------------------------------------------------------------

# The physical argument of an equation with a closure unknown, under a key that
# no state function reads. Added after the arguments of all equations are
# merged, so it never leaks to a sibling equation.
@inline _with_secant_argument(args, ::Any, x) = args
@inline _with_secant_argument(args, closure::Int, x) = merge(args, (; secant_argument = x[closure]))

"""
    _evaluate_leaf(fn, r, args)

`fn(r, args)`, unless `args` carries a `secant_argument` and `r` uses a secant
viscosity. Then `r` is evaluated at that physical argument `s` and its
contribution rescaled onto the equation's own (shifted) unknown `v`:
`fn(r; k = |s|) * v / |s|`, where `k` is the key of the unknown.
"""
@inline _evaluate_leaf(fn::F, r, args) where {F} =
    _evaluate_leaf(Val(_uses_secant(r) && hasfield(typeof(args), :secant_argument)), fn, r, args)

@inline _evaluate_leaf(::Val{false}, fn::F, r, args) where {F} = fn(r, args)

@inline function _evaluate_leaf(::Val{true}, fn::F, r, args) where {F}
    s = abs(args.secant_argument)
    k = first(keys(differentiable_kwargs(fn)))
    v = getfield(args, k)
    return fn(r, merge(args, NamedTuple{(k,)}((s,)))) * v * safe_inv(s)
end

# -----------------------------------------------------------------------
# Closure residuals
# -----------------------------------------------------------------------

"""
    add_closure_residuals(residual, eqs, x, E, others)

Add to `residual` the residuals of the closure equations, which compare each
closure unknown with the invariant of the physical tensor it stands for. `E` is
the assembled strain rate `ε + H`.
"""
@inline add_closure_residuals(residual::NTuple{N, Any}, eqs::NTuple{N, CompositeEquation}, x, E, others) where {N} =
    maptuple((r, eq) -> r + closure_residual(eq, x, E, others), residual, eqs)

@inline closure_residual(::CompositeEquation, x, E, others) = false

@inline function closure_residual(eq::CompositeEquation{false, T, F, R, RT, <:BlockClosure}, x, E, others) where {T, F, R, RT}
    idx = eq.closure
    seed = E .* zero(eltype(x))
    args = merge((; ε = x[idx.block]), others)
    phys = (; e = _x_at(x, idx.e), t = maptuple(i -> _x_at(x, i), idx.t), el = eq.el_number)
    block = eq.rheology
    state = _block_state(block.leafs, block.branches, _backstresses(others, seed), args, idx.τ0_start, phys, seed)
    return _closure_residual(eq.fn, idx, eq.self, state, _direction(E), x)
end

@inline function _closure_residual(::typeof(compute_block_strain_rate_invariant), idx, self, state, n, x)
    η_KV, W = state
    return x[self] - differentiable_second_invariant(x[idx.block] .* n .- W ./ (2 * η_KV))
end

@inline function _closure_residual(::typeof(compute_subbranch_stress_invariant), idx::BlockClosure{J}, self, state, n, x) where {J}
    η_KV, W, ηM, S = state
    return x[self] - differentiable_second_invariant(x[idx.subs[J]] .* n .- (ηM[J] / η_KV) .* W .+ S[J])
end

# Unit direction of the assembled strain rate; a signed scalar strain rate has
# none, its unknowns carry the sign.
@inline _direction(E::Number) = one(E)
@inline _direction(E::NTuple{1}) = (one(E[1]),)
@inline _direction(E::NTuple) = E .* safe_inv(differentiable_second_invariant(E))

# `others.τ0`, or zero backstresses of the shape of `seed` when `others` has none.
struct ZeroBackstresses{T}
    value::T
end
Base.getindex(z::ZeroBackstresses, ::Int) = z.value

@inline _backstresses(others, seed) = hasfield(typeof(others), :τ0) ? others.τ0 : ZeroBackstresses(seed)
