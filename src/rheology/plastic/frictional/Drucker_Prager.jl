"""
    DruckerPrager{T} <: AbstractPlasticity

Represents the Drucker-Prager plasticity model for pressure-dependent yielding.

# Fields
- `C::T`: The cohesion parameter.
- `ϕ::T`: The friction angle (in degrees).
- `ψ::T`: The dilatancy angle (in degrees).
- `η_vp::T`: The Duvaut-Lions viscoplastic regularisation viscosity.
"""
struct DruckerPrager{T} <: AbstractPlasticity
    C::T
    ϕ::T # in degrees for now
    ψ::T # in degrees for now
    η_vp::T # regularisation viscosity
    sinϕ::T
    sinψ::T
    cosϕ::T
    cosψ::T

    function DruckerPrager(C::T, ϕ::T, ψ::T, η_vp::T) where {T}
        sinϕ, cosϕ = sincosd(ϕ)
        sinψ, cosψ = sincosd(ψ)
        return new{T}(C, ϕ, ψ, η_vp, sinϕ, sinψ, cosϕ, cosψ)
    end
end
function DruckerPrager(C, ϕ, ψ)
    C, ϕ, ψ = promote(C, ϕ, ψ)
    return DruckerPrager(C, ϕ, ψ, one(C))
end
DruckerPrager(args::Vararg{Any, 4}) = DruckerPrager(promote(args...)...)

@inline _isvolumetric(::DruckerPrager) = false

@inline series_state_functions(::DruckerPrager) = (compute_strain_rate, compute_lambda)

@inline parallel_state_functions(::DruckerPrager) = (compute_stress, compute_pressure, compute_lambda_parallel, compute_plastic_strain_rate)

@inline function compute_strain_rate(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r::DruckerPrager; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    return ε_pl / 2
end

@inline function compute_volumetric_strain_rate(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    θ_pl = compute_volumetric_plastic_strain_rate(r::DruckerPrager; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    return -θ_pl
end

@inline function compute_lambda(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    return F - λ * r.η_vp
end

@inline function compute_lambda_parallel(r::DruckerPrager; τ_pl = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ_pl, P)
    return F - λ * r.η_vp
end

# special plastic helper functions
function compute_F(r::DruckerPrager, τ, P)
    F = (τ - P * r.sinϕ - r.C * r.cosϕ)
    return F * (F > -1.0e-8)
end
compute_Q(r::DruckerPrager, τ, P) = τ - P * r.sinψ

@inline function compute_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ - ε
end

@inline function compute_volumetric_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return -λ * r.sinψ - θ
end

@inline compute_plastic_stress(r::DruckerPrager; τ_pl = 0, kwargs...) = τ_pl
@inline compute_stress(r::DruckerPrager; τ_pl = 0, kwargs...) = τ_pl

@inline compute_viscosity(r::DruckerPrager; kwargs...) = Inf
