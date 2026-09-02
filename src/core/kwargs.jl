"""
    history_kwargs(r)

Return the names of `others` fields that are interpreted as element-local
history tuples for rheology `r`.
"""
history_kwargs(::AbstractElasticity) = (:τ0, :P0)
history_kwargs(::AbstractViscosity) = (:d,)
history_kwargs(::AbstractPlasticity) = ()

"""
    differentiable_kwargs(T, fn)
    differentiable_kwargs(fn)

Return a zero-valued `NamedTuple` describing the unknown solved by a state
function. These keys define how entries of `x` are exposed to rheology methods.
"""
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_strain_rate)) where {T} = (; τ = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_volumetric_strain_rate)) where {T} = (; P = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_lambda)) where {T} = (; λ = zero(T)) # τ = zero(T), P = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_lambda_parallel)) where {T} = (; λ = zero(T)) # τ = zero(T), P = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_stress)) where {T} = (; ε = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_pressure)) where {T} = (; θ = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_plastic_strain_rate)) where {T} = (; τ_pl = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_plastic_stress)) where {T} = (; τ_pl = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_volumetric_plastic_strain_rate)) where {T} = (; τ_pl = zero(T), P_pl = zero(T))
@inline differentiable_kwargs(fun::F) where {F <: Function} = differentiable_kwargs(Float64, fun)

@inline differentiable_kwargs(::Tuple{}) = (;)
@inline differentiable_kwargs(funs::NTuple{N, Any}) where {N} = differentiable_kwargs(Float64, funs)

@generated function differentiable_kwargs(::Type{T}, funs::NTuple{N, Any}) where {N, T}
    return quote
        @inline
        Base.@nexprs $N i -> nt_i = differentiable_kwargs($T, funs[i])
        Base.@ncall $N merge nt
    end
end
differentiable_kwargs(::Type{T}, funs::NTuple{0, Any}) where {T} = (;)
