# State functions shared by every AbstractCapPlasticity element. A subtype
# supplies only its struct, compute_F, and compute_Q; everything below is
# expressed in terms of those two and so needs no per-model definition.

@inline _isvolumetric(::AbstractCapPlasticity) = true

@inline series_state_functions(::AbstractCapPlasticity) = (compute_strain_rate, compute_lambda, compute_volumetric_strain_rate)
@inline parallel_state_functions(::AbstractCapPlasticity) = compute_stress, compute_pressure, compute_lambda, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate

@inline function compute_strain_rate(r::AbstractCapPlasticity; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    F = compute_F(r, τ, P)
    return ε_pl / 2 * (F > -1.0e-8)
end
@inline function compute_volumetric_strain_rate(r::AbstractCapPlasticity; τ = 0, λ = 0, P = 0, kwargs...)
    θ_pl = compute_volumetric_plastic_strain_rate(r; τ_pl = τ, λ = λ, P_pl = P, θ = 0, kwargs...)
    F = compute_F(r, τ, P)
    return θ_pl * (F > -1.0e-8)
end

# Consistency residual: `F = λ·η_vp` above yield, `λ·η_vp = 0` below it.
@inline function compute_lambda(r::AbstractCapPlasticity; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    yielding = F > -1.0e-8
    return -F * yielding + λ * _lambda_below_yield(r, yielding)
end

# Below yield an unregularized cap (`η_vp = 0`) would leave `λ` undetermined,
# so the factor falls back to `oneunit` there to pin it at zero.
@inline _lambda_below_yield(r::AbstractCapPlasticity, yielding) =
    (iszero(r.η_vp) && !yielding) ? oneunit(r.η_vp) : r.η_vp

@inline compute_stress(r::AbstractCapPlasticity; τ_pl = 0, kwargs...) = τ_pl
@inline compute_pressure(r::AbstractCapPlasticity; P_pl = 0, kwargs...) = P_pl

# The flow rule is the gradient of the potential Q, taken by automatic
# differentiation so a subtype needs only to define Q itself.
@inline function compute_plastic_strain_rate(r::AbstractCapPlasticity; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ * ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) - 0 * ε
end

@inline function compute_volumetric_plastic_strain_rate(r::AbstractCapPlasticity; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return -λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) - 0 * θ
end
