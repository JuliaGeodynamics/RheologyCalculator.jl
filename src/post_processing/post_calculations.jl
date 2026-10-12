# Postprocessing calculations (not needed as part of local iterations)

"""
    compute_stress_elastic(c, xnew, others)
    compute_stress_elastic(eqs, xnew, others)

Return a flat tuple with the updated stresses of all elastic elements present in
`c` or equation tuple `eqs`.

This is a post-processing helper for advancing elastic history fields after a
successful solve. For equations written in strain-rate form, the corresponding
entry of `xnew` is returned; for equations written in stress form, the elastic
state function is evaluated from the local arguments.
"""
compute_stress_elastic(c::AbstractCompositeModel, xnew, others) = compute_stress_elastic(generate_equations(c), xnew, others)

@inline isglobal(::CompositeEquation{B}) where {B} = Val(B)

@inline function compute_stress_elastic(eqs::NTuple{N, CompositeEquation}, xnew, others) where {N}
    args_all = generate_args_template(eqs, xnew, others)
    τ_elastic = maptuple(eqs, args_all) do eq, args
        (; self, rheology, fn, el_number) = eq
        _compute_stress_elastic(rheology, self, fn, el_number, xnew, others, args, isglobal(eq))
    end
    return superflatten(τ_elastic)
end

@inline _compute_stress_elastic(rheology::NTuple{N, Any}, self, fn::F, el_number, xnew, others, args, eqglobal) where {N, F} =
    superflatten(maptuple((r, number) -> _compute_stress_elastic1(r, self, fn, number, xnew, others, args, eqglobal), rheology, el_number))

# A closure equation holds a whole block and no elastic stress of its own.
@inline _compute_stress_elastic(::AbstractCompositeModel, ::Vararg{Any, N}) where {N} = ()

@inline _compute_stress_elastic1(::T, ::Vararg{Any, N}) where {T, N} = ()
function _compute_stress_elastic1(r::AbstractElasticity, self, fn::F, number, xnew, others, args, eqglobal) where {F}

    @inline f(::F, ::Vararg{Any, N}) where {F, N} = ()

    # Direct elastic leaf of the outer SeriesModel: the effective strain-rate
    # correction (see effective_strain_rate_correction / docs/derivations/tensor_reduction.typ)
    # already absorbed τ0 into the global solve, so x[self] IS the physical
    # spring stress.
    @inline f(::typeof(compute_strain_rate), r, others, args, number, self, xnew, ::Val{true}) = xnew[self]

    # Elastic element of a sub-branch: x[self] is the shifted sub-branch stress,
    # not the physical spring stress. Springs inside blocks are recovered by
    # `_spring_stresses` instead.
    @inline function f(::typeof(compute_strain_rate), r, others, args, number, self, xnew, ::Val{false})
        return ()
    end

    # compute_stress is never IsGlobal: the outer composite is always a
    # SeriesModel, whose global equations are compute_strain_rate /
    # compute_volumetric_strain_rate; compute_stress equations only arise
    # inside ParallelModel branches, which are never IsGlobal = true. Falls
    # through to the generic fallback above for both Val{true} and Val{false}.
    @inline function f(::typeof(compute_stress), r, others, args, number, self, xnew, ::Val{false})
        return ()
    end

    return f(fn, r, others, args, number, self, xnew, eqglobal)
end

"""
    compute_stress_elastic(c::SeriesModel, x, others)
    compute_stress_elastic(c::SeriesModel, x, ε, others)

Return the stress of every elastic element of `c` after a solve, one entry per
spring in the order of `τ0`: the second invariant of the spring stress tensor
of [`elastic_stress_history`](@ref), signed for a scalar strain rate.

A spring inside a `ParallelModel` block is not coaxial with the outer stress
once the block carries a tensor history, so recovering it needs the applied
strain rate `ε`. The form without `ε` therefore applies only when every spring
is a leaf of the outer `SeriesModel`, or when `others.τ0` holds scalars;
otherwise it throws.
"""
compute_stress_elastic(c::SeriesModel, x::AbstractVector, others) =
    _compute_stress_elastic_without_ε(Val(_has_block_springs(c)), c, x, others)

compute_stress_elastic(c::SeriesModel, x::AbstractVector, ε, others) =
    maptuple(second_invariant_value, _spring_stress_tensors(c, _static_solution(c, x), ε, others))

_compute_stress_elastic_without_ε(::Val{false}, c, x, others) = compute_stress_elastic(generate_equations(c), x, others)

function _compute_stress_elastic_without_ε(::Val{true}, c, x, others)
    _scalar_history(others) || _throw_stress_elastic_needs_ε()
    # A scalar history has no direction: the spring stresses are signed scalars.
    xs = _static_solution(c, x)
    return _spring_stresses(c, xs, one(eltype(xs)), others)
end

@inline _has_block_springs(c::SeriesModel) = foldtuple(|, false, _iselastic, c.branches)

@inline _scalar_history(others) = hasfield(typeof(others), :τ0) ? _scalar_history_entries(others.τ0) : true
@inline _scalar_history_entries(::Tuple{Vararg{Number}}) = true
@inline _scalar_history_entries(::Tuple) = false

@noinline function _throw_stress_elastic_needs_ε()
    throw(
        ArgumentError(
            "the springs inside a ParallelModel block carry a tensor history, so their stresses " *
                "depend on the applied strain rate; use `compute_stress_elastic(c, x, ε, others)` " *
                "or `elastic_stress_history(c, x, ε, τ0, others)`"
        )
    )
end

# The solver vector as an `SVector`, as the residual is evaluated on.
@inline function _static_solution(c, x::AbstractVector)
    n = length(generate_equations(c))
    length(x) == n || throw(DimensionMismatch("x has $(length(x)) entries but the solver vector of the composite has $n"))
    return SVector{n}(ntuple(i -> x[i], n))
end
@inline _static_solution(c, x::SVector) = x

"""
    elastic_stress_history(c::SeriesModel, x, ε, τ0, others)

Advance the elastic backstress history of composite `c` by one step: return the
stress tensor of every spring, in the order of `τ0`, for use as the next step's
`τ0`. `x` is the converged solver vector of the step solved with strain rate `ε`
and history `τ0`; each entry has the shape of `ε`, a Voigt tuple or a signed
scalar.

With `E = ε + H` the assembled strain rate and `n = E / E_II` its direction
(`n = 1` for a scalar `ε`), the outer stress is `τ = τ_II n`, and for each
`ParallelModel` block `b` with shifted strain-rate unknown `x̃_b`

```
ε̇_b = x̃_b n - W_b / (2 η_KV,b)
```

is its physical strain rate. The springs carry

- `τ`, for a spring of the outer `SeriesModel`;
- `τ0_s + 2 G_s Δt ε̇_b`, for a direct spring of block `b`;
- `τ_j = 2 η_M,j ε̇_b + S_j`, for every spring of sub-branch `j` of block `b`,

with `η_KV,b`, `W_b`, `η_M,j` and `S_j` the coefficients of the block, evaluated
as in the solve (see [`effective_strain_rate_correction`](@ref)).
"""
elastic_stress_history(c::SeriesModel, x::AbstractVector, ε, τ0, others) =
    _spring_stress_tensors(c, _static_solution(c, x), ε, merge(others, (; τ0)))

@inline function _spring_stress_tensors(c::SeriesModel, x::SVector, ε, others)
    eqs = generate_equations(c)
    n = _direction(assembled_strain_rate(c, eqs, x, ε, others))
    return _spring_stresses(c, x, n, others)
end

# The spring stresses for outer-stress direction `n`: a unit tensor, or `1` for a
# signed scalar.
@inline function _spring_stresses(c::SeriesModel, x::SVector, n, others)
    eqs = generate_equations(c)
    τ = x[_global_deviatoric_index(eqs)] .* n
    seed = n .* zero(eltype(x))
    outer = ntuple(_ -> τ, Val(count_elastic(c.leafs)))
    εb = _block_strain_rates(eqs, x, c.branches)
    phys = _block_closure_values(eqs, x, c.branches, global_eltype_numbering(c)[2])
    blocks = _block_spring_stresses(c.branches, _backstresses(others, seed), others, count_elastic(c.leafs), εb, phys, n, seed)
    return (outer..., blocks...)
end

@generated function _global_deviatoric_index(eqs::NTuple{N, Any}) where {N}
    k = findfirst(E -> E.parameters[1] === true && E.parameters[3] === typeof(compute_strain_rate), collect(eqs.parameters))
    k === nothing && return :(error("composite has no global deviatoric equation"))
    return :($k)
end

# One tuple of spring stresses per branch, concatenated in branch order; each
# branch reads its backstresses from a literal τ0 offset.
@generated function _block_spring_stresses(branches::NTuple{Nb, Any}, τ0, others, offset, εb, phys, n, seed) where {Nb}
    offsets = _branch_tau0_offsets(branches)
    parts = [
        :(_branch_spring_stresses(iselastic(branches[$i]), branches[$i], τ0, others, offset + $(offsets[i]), εb[$i], phys[$i], n, seed)...)
            for i in 1:Nb
    ]
    return quote
        @inline
        $(Expr(:tuple, parts...))
    end
end

@inline _branch_spring_stresses(::Val{false}, branch, τ0, others, el_idx_start, εb, phys, n, seed) = ()

@inline function _branch_spring_stresses(::Val{true}, branch::ParallelModel, τ0, others, el_idx_start, εb, phys, n, seed)
    args, ε_block, ηM, S = _physical_block_state(branch, τ0, others, el_idx_start, εb, phys, n, seed)
    return _block_springs(branch.leafs, branch.branches, τ0, args, el_idx_start, ε_block, ηM, S)
end

# The physical strain rate `ε_block = x̃_b n - W / (2 η_KV)` of a block that
# carries elastic history, with the coefficients it was computed from.
@inline function _physical_block_state(branch::ParallelModel, τ0, others, el_idx_start, εb, phys, n, seed)
    args = merge((; ε = εb), others)
    η_KV, W, ηM, S = _block_state(branch.leafs, branch.branches, τ0, args, el_idx_start, phys, seed)
    ε_block = εb .* n .- W ./ (2 * η_KV)
    return args, ε_block, ηM, S
end

"""
    physical_solution(c::SeriesModel, x, ε, others)

The solver vector `x` with the shifted unknowns of every `ParallelModel` block
that carries elastic history replaced by physical invariants: the block strain
rate by `‖ε̇_b‖` and each sub-branch stress by `‖τ_j‖` (signed values for a
scalar `ε`), with the tensors of [`elastic_stress_history`](@ref). Other
entries are unchanged. Element state functions evaluated on this vector give
each element's physical strain rate and stress invariants.
"""
function physical_solution(c::SeriesModel, x::SVector, ε, others)
    eqs = generate_equations(c)
    n = _direction(assembled_strain_rate(c, eqs, x, ε, others))
    seed = n .* zero(eltype(x))
    iblock = block_equation_indices(eqs, c.branches)
    phys = _block_closure_values(eqs, x, c.branches, global_eltype_numbering(c)[2])
    return _physical_blocks(c.branches, x, iblock, _backstresses(others, seed), others, count_elastic(c.leafs), phys, n, seed, global_eltype_numbering(c)[2])
end

@generated function _physical_blocks(branches::NTuple{Nb, Any}, x, iblock, τ0, others, offset, phys, n, seed, el_nums) where {Nb}
    offsets = _branch_tau0_offsets(branches)
    stmts = [
        :(x = _physical_block(iselastic(branches[$i]), branches[$i], x, iblock[$i], τ0, others, offset + $(offsets[i]), phys[$i], n, seed, el_nums[$i]))
            for i in 1:Nb
    ]
    return quote
        @inline
        $(stmts...)
        return x
    end
end

@inline _physical_block(::Val{false}, branch, x, ib, τ0, others, el_idx_start, phys, n, seed, el_num) = x

@inline function _physical_block(::Val{true}, branch::ParallelModel, x, ib, τ0, others, el_idx_start, phys, n, seed, el_num)
    _, ε_block, ηM, S = _physical_block_state(branch, τ0, others, el_idx_start, x[ib], phys, n, seed)
    isubs = branch_equation_indices(branch, compute_stress, ib - 1, el_num)
    x = setindex(x, differentiable_second_invariant(ε_block), ib)
    return foldtuple((xi, (j, i)) -> setindex(xi, differentiable_second_invariant((2 * ηM[j]) .* ε_block .+ S[j]), i), x, identity, tuple_pairs(isubs))
end

# `((1, i₁), (2, i₂), …)` for `(i₁, i₂, …)`.
@inline tuple_pairs(t::NTuple{N, Any}) where {N} = ntuple(j -> (j, t[j]), Val(N))

# Direct springs first, then the springs of each sub-branch, matching the τ0
# order of `_block_state`.
@generated function _block_springs(
        leafs::NTuple{N, AbstractRheology}, subs::NTuple{Ns, Any}, τ0, args, el_idx_start, ε_block, ηM, S
    ) where {N, Ns}
    direct, per_sub = _elastic_source_positions(leafs, subs)
    entries = Any[]
    el_count = 0
    for pos in direct
        el_count += 1
        push!(entries, :(τ0[el_idx_start + $el_count] .+ (2 * compute_viscosity(leafs[$pos], args)) .* ε_block))
    end
    for j in 1:Ns
        isempty(per_sub[j]) && continue
        τ_j = Symbol(:τ_, j)
        for _ in per_sub[j]
            push!(entries, τ_j)
        end
    end
    sub_stresses = [:($(Symbol(:τ_, j)) = (2 * ηM[$j]) .* ε_block .+ S[$j]) for j in 1:Ns if !isempty(per_sub[j])]
    return quote
        @inline
        $(sub_stresses...)
        $(Expr(:tuple, entries...))
    end
end

"""
    compute_pressure_elastic(c, xnew, others)
    compute_pressure_elastic(eqs, xnew, others)

Return a flat tuple with the updated pressures of all elastic elements present
in `c` or equation tuple `eqs`.

This is the volumetric counterpart of `compute_stress_elastic` and is intended
for updating pressure history after a successful solve.
"""
@inline function compute_pressure_elastic(eqs::NTuple{N, CompositeEquation}, xnew, others) where {N}
    args_all = generate_args_template(eqs, xnew, others)
    P_elastic = maptuple(eqs, args_all) do eq, args
        (; self, rheology, fn, el_number) = eq
        _compute_pressure_elastic(rheology, self, fn, el_number, xnew, others, args)
    end
    return superflatten(P_elastic)
end
compute_pressure_elastic(c::AbstractCompositeModel, xnew, others) = compute_pressure_elastic(generate_equations(c), xnew, others)

@inline _compute_pressure_elastic(rheology::NTuple{N, Any}, self, fn::F, el_number, xnew, others, args) where {N, F} =
    superflatten(maptuple((r, number) -> _compute_pressure_elastic1(r, self, fn, number, xnew, others, args), rheology, el_number))

@inline _compute_pressure_elastic1(::T, ::Vararg{Any, N}) where {T, N} = ()
function _compute_pressure_elastic1(r::AbstractElasticity, self, fn::F, number, xnew, others, args) where {F}

    @inline f(::F, ::Vararg{Any, N}) where {F, N} = ()
    @inline f(::typeof(compute_volumetric_strain_rate), r, others, args, number, self, xnew) = xnew[self]

    @inline function f(::typeof(compute_pressure), r, others, args, number, self, xnew)
        # compute elastic pressure from local volumetric strain rate
        keys_hist = history_kwargs(r)
        args_local = extract_local_kwargs(others, keys_hist, number)
        args_combined = merge(args, args_local)
        return compute_pressure(r, args_combined)
    end

    return f(fn, r, others, args, number, self, xnew)
end

"""
    volumetric_plastic_strain_rate(c::AbstractCompositeModel, x, others)

Return the volumetric plastic strain rate `θ_pl` that composite `c` produces at
solver vector `x`, which is the plastic contribution a Stokes solver adds to
mass conservation.

Only dilatant elements contribute, and their contribution is a function of the
solved plastic multiplier, so it cannot be recovered by calling
[`compute_volumetric_strain_rate`](@ref) on the composite with a prescribed
pressure: that call leaves `λ` at its zero default and returns no plastic part.

`x` and `others` must be the solver vector and the auxiliary fields passed to
[`solve`](@ref).

# Example
```julia
c = SeriesModel(LinearViscosity(2.5e3), Elasticity(1e3, 1e3), DruckerPrager(2.0, 30.0, 30.0))
others = (; dt = 0.1, τ0 = ((0.0, 0.0, 0.0),), P0 = (0.0,))
sol = solve(c, initial_guess_x(c, (; ε = 1e-2, θ = 0.0), (; τ = 0.0), others), (; ε = (1e-2, -1e-2, 0.0), θ = 0.0), others)
volumetric_plastic_strain_rate(c, sol, others)
```
"""
function volumetric_plastic_strain_rate(c::AbstractCompositeModel, sol::AbstractVector, others)
    eqs = generate_equations(c)
    n = length(eqs)
    x = SVector{n}(ntuple(i -> sol[i], n))
    args = generate_args_template(eqs, x, others)
    return _volumetric_plastic_strain_rate(eqs, args, others)
end

@generated function _volumetric_plastic_strain_rate(eqs::NTuple{N, CompositeEquation}, args, others) where {N}
    N == 0 && return :(0.0)
    return quote
        @inline
        sum(Base.@ntuple $N i -> _volumetric_plastic_equation(eqs[i].fn, eqs[i], args[i], others))
    end
end

@inline _volumetric_plastic_equation(::F, _eq, _args, _others) where {F} = 0.0

# A plastic element owns a volumetric equation in exactly one of these two
# forms, series or parallel, so summing over both counts it once.
@inline function _volumetric_plastic_equation(
        fn::Union{typeof(compute_volumetric_strain_rate), typeof(compute_volumetric_plastic_strain_rate)},
        eq, args, others,
    )
    values = evaluate_state_function_perleaf(fn, eq.rheology, args, others, eq.el_number)
    return _sum_plastic(eq.rheology, values)
end

@inline _sum_plastic(::Tuple{}, ::Tuple{}) = 0.0

@generated function _sum_plastic(rheology::NTuple{N, AbstractRheology}, values) where {N}
    return quote
        @inline
        sum(Base.@ntuple $N i -> _plastic_contribution(rheology_category(rheology[i]), values[i]))
    end
end

@inline _plastic_contribution(::Val{:plastic}, v) = v
@inline _plastic_contribution(::Val, v) = zero(v)
