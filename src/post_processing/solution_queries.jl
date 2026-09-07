# Semantic access to entries of the local solution vector.
#
# `x_keys` labels each unknown by its *kind* -- `:τ`, `:ε`, `:P`, `:λ` -- and not by
# where it sits in the network. A nested composite therefore repeats keys: a Maxwell
# element inside a parallel branch yields `(:τ, :ε, :τ)`, a two-unit Kelvin chain
# yields `(:τ, :ε, :ε)`. Asking for "the" stress by symbol is ambiguous on any such
# model, so the primary quantities are selected by *topology* instead: the unique
# global equation of the appropriate kind.
#
# All queries take the composite model directly, matching the `f(c, x, ...)` grammar
# of `solve`, `compute_residual` and `tangent`. Solution arguments are annotated
# `::AbstractVector` rather than `::SVector` so a wrapper type around the solution
# vector dispatches to the same method.

"""
    solution_values(c::AbstractCompositeModel, x::AbstractVector, ::Val{k})

Return every entry of the solution vector `x` whose key is `k`, as a tuple, in
solver order. Repeated keys are preserved; a key the model does not use returns
an empty tuple.

Use this when a model may carry several unknowns of the same kind and you want
all of them. For the single global stress or pressure of a model, prefer
[`primary_deviatoric_stress`](@ref) or [`primary_pressure`](@ref), which select
by topology rather than by symbol.

```julia
c = SeriesModel(v, ParallelModel(SeriesModel(v2, e), v3))   # x_keys -> (:τ, :ε, :τ)
solution_values(c, x, Val(:τ))                              # -> (x[1], x[3])
```
"""
@inline solution_values(c::AbstractCompositeModel, x::AbstractVector, ::Val{k}) where {k} =
    _values_with_key(x, Val(x_keys(c)), Val(k))

# `ks` is lifted into a type parameter so the matching positions are resolved at
# compile time and the body reduces to a literal tuple construction.
@generated function _values_with_key(x::AbstractVector, ::Val{ks}, ::Val{k}) where {ks, k}
    idx = [i for i in eachindex(ks) if ks[i] === k]
    return quote
        @inline
        @inbounds ($([:(x[$i]) for i in idx]...),)
    end
end

@inline _isglobal(::CompositeEquation{B}) where {B} = B

# Both helpers resolve entirely from the equation *types*: `CompositeEquation`
# carries its globality as its first type parameter and its state function as its
# third, so nothing below survives into the generated code as a runtime test.
@generated function _count_global(eqs::NTuple{N, CompositeEquation}, ::F) where {N, F}
    n = count(i -> eqs.parameters[i].parameters[1] === true &&
            eqs.parameters[i].parameters[3] === F, 1:N)
    return :($n)
end

@generated function _first_global(eqs::NTuple{N, CompositeEquation}, ::F) where {N, F}
    i = findfirst(j -> eqs.parameters[j].parameters[1] === true &&
            eqs.parameters[j].parameters[3] === F, 1:N)
    return i === nothing ? :(0) : :($i)
end

@inline function _unique_global_index(eqs, fn::F, what) where {F}
    n = _count_global(eqs, fn)
    isone(n) || throw(
        ArgumentError(
            "composite has $n global $what equations; exactly one is required. " *
                "Use `solution_values` to retrieve every entry of a repeated key, " *
                "or supply a model-specific method."
        )
    )
    return _first_global(eqs, fn)
end

"""
    primary_stress_index(c::AbstractCompositeModel)

Position in the solution vector of the model's global deviatoric stress: the
unique equation that is global and whose state function is `compute_strain_rate`.

Throws an `ArgumentError` unless there is exactly one such equation. Selecting on
globality alone is not enough -- a compressible model has both a deviatoric and a
volumetric global equation -- so the state function is part of the test.
"""
@inline primary_stress_index(c::AbstractCompositeModel) =
    _unique_global_index(generate_equations(c), compute_strain_rate, "deviatoric stress")

"""
    primary_pressure_index(c::AbstractCompositeModel)

Position of the model's global volumetric pressure unknown, or `nothing` when the
model has no volumetric equation. Throws if several are found.
"""
@inline function primary_pressure_index(c::AbstractCompositeModel)
    eqs = generate_equations(c)
    iszero(_count_global(eqs, compute_volumetric_strain_rate)) && return nothing
    return _unique_global_index(eqs, compute_volumetric_strain_rate, "volumetric pressure")
end

"""
    primary_deviatoric_stress(c::AbstractCompositeModel, x::AbstractVector)

The model's global deviatoric stress invariant: the value of the unique global
`compute_strain_rate` equation.

This is the quantity a momentum equation consumes. It is *not* simply the first
`:τ` in `x_keys` -- a nested composite has several `:τ` unknowns, one per
`SeriesModel` node, and only one of them is the stress carried by the composite
as a whole. Throws an `ArgumentError` when the model has zero or several global
deviatoric equations.
"""
@inline primary_deviatoric_stress(c::AbstractCompositeModel, x::AbstractVector) =
    @inbounds x[primary_stress_index(c)]

"""
    primary_pressure(c::AbstractCompositeModel, x::AbstractVector; fallback = nothing)

The model's global volumetric pressure unknown, or `fallback` when the model is
incompressible and carries no volumetric equation.

`fallback` is a keyword because it is caller policy rather than solution data: an
incompressible composite has no pressure to report, and what should stand in its
place -- the trial pressure, zero, or nothing at all -- is the caller's decision.
"""
@inline function primary_pressure(c::AbstractCompositeModel, x::AbstractVector; fallback = nothing)
    i = primary_pressure_index(c)
    i === nothing && return fallback
    return @inbounds x[i]
end

"""
    plastic_multipliers(c::AbstractCompositeModel, x::AbstractVector)

Every plastic multiplier in the solution vector, as a tuple in solver order.
Returns an empty tuple for a model without plasticity, rather than throwing:
absence of plasticity is an ordinary case, not an error.
"""
@inline plastic_multipliers(c::AbstractCompositeModel, x::AbstractVector) =
    solution_values(c, x, Val(:λ))
