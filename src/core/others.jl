"""
    superflatten(x)

Recursively flatten nested tuples and return a single flat tuple. Non-tuple
values are wrapped as one-element tuples.
"""
@inline superflatten(t::NTuple{N, Any}) where {N} = superflatten(first(t))..., superflatten(Base.tail(t))...
@inline superflatten(::Tuple{}) = ()
@inline superflatten(x) = (x,)

"""
    isvolumetric(c)

Return `Val(true)` if a rheology or composite contains any volumetric state
functions, otherwise `Val(false)`.
"""
isvolumetric(c::AbstractCompositeModel) = Val(_isvolumetric(c))
isvolumetric(c::AbstractRheology) = Val(_isvolumetric(c))

@inline _isvolumetric(r::NTuple{N, AbstractRheology}) where {N} = foldtuple(|, false, _isvolumetric, r)

@inline _isvolumetric(::AbstractRheology) = false
# @inline _isvolumetric(::Elasticity) = true
# @inline _isvolumetric(::BulkElasticity) = true
# @inline _isvolumetric(::BulkViscosity) = true
# @inline _isvolumetric(c::AbstractCompositeModel) = _isvolumetric(c.leafs)
@inline _isvolumetric(::Tuple{}) = false

_isvolumetric(c::AbstractCompositeModel) = _isvolumetric(c.leafs, c.branches)

@inline _isvolumetric(leafs, branches::Tuple) =
    _isvolumetric(leafs) | foldtuple(|, false, _isvolumetric, branches)

# Sparsity-detection tracers (which carry no primal value)
# can override those functions in the SparseConnectivityTracer extension. The guards only
# ever protected the value against Inf/NaN. The dependency pattern is the
# same with or without them. Float64 behaviour is unchanged.
@inline safe_inv(v) = iszero(primal(v)) ? zero(v) : inv(v)
@inline safe_inv_one(v) = iszero(primal(v)) ? one(v) : inv(v)

# Value of `x` with every layer of dual partials stripped.
# `iszero` on a `ForwardDiff.Dual` also tests the partials, so a zero primal
# carrying a nonzero derivative (a strain rate of zero seeded for differentiation)
# would otherwise take the `inv` branch and produce `Inf` and `NaN` partials.
@inline primal(x) = x
@inline primal(x::ForwardDiff.Dual) = primal(ForwardDiff.value(x))
