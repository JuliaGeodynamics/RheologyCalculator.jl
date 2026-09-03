import ..RheologyCalculator: series_state_functions, parallel_state_functions, _isvolumetric
import ..RheologyCalculator: compute_strain_rate, compute_stress, compute_pressure, compute_volumetric_strain_rate
import ..RheologyCalculator: compute_plastic_strain_rate, compute_plastic_stress, compute_volumetric_plastic_strain_rate
import ..RheologyCalculator: compute_lambda, compute_lambda_parallel
import ..RheologyCalculator: compute_viscosity, compute_viscosity_series, compute_viscosity_parallel

"""
    DruckerPrager{T} <: AbstractPlasticity

Represents the Drucker-Prager plasticity model for pressure-dependent yielding.

# Fields
- `C::T`: The cohesion parameter.
- `ϕ::T`: The friction angle (in degrees).
- `ψ::T`: The dilatancy angle (in degrees).
"""
struct DruckerPrager{T} <: AbstractPlasticity
    C::T
    ϕ::T # in degrees for now
    ψ::T # in degrees for now
    sinϕ::T
    sinψ::T
    cosϕ::T
    cosψ::T

    function DruckerPrager(C::T, ϕ::T, ψ::T) where T
        sinϕ, cosϕ = sincosd(ϕ)
        sinψ, cosψ = sincosd(ψ)
        new{T}(C, ϕ, ψ, sinϕ, sinψ, cosϕ, cosψ)
    end
end
DruckerPrager(args::Vararg{Any, 3}) = DruckerPrager(promote(args...)...)

@inline _isvolumetric(::DruckerPrager) = false

@inline function series_state_functions(r::DruckerPrager)
    # we need to check whether this allocates
    # if r.ψ == 0
        return (compute_strain_rate, compute_lambda)
    # else
        # return (compute_strain_rate, compute_volumetric_strain_rate, compute_lambda)
    # end
end

@inline function parallel_state_functions(r::DruckerPrager)
    # if r.ψ == 0
        return (compute_stress, compute_pressure, compute_lambda_parallel, compute_plastic_strain_rate)
    # else
        # return (compute_stress, compute_pressure, compute_lambda_parallel, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate)
    # end
end

@inline function compute_strain_rate(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r::DruckerPrager; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    F = compute_F(r, τ, P)
    return ε_pl/2#*(F > -1e-8)
end

@inline function compute_volumetric_strain_rate(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    θ_pl = compute_volumetric_plastic_strain_rate(r::DruckerPrager; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    F = compute_F(r, τ, P)
    return -θ_pl#*(F > -1e-8)
end

@inline function compute_lambda(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    η_χ = 1.0  # Lagrange multiplier, value doesn't matter
    # return F*(F>-1e-8) - λ * η_χ #* (F < -1e-8)
    return F - λ * η_χ #* (F < -1e-8)
end

@inline function compute_lambda_parallel(r::DruckerPrager; τ_pl = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ_pl, P)
    η_χ = 1.0  # Lagrange multiplier, value doesn't matter
    return F - λ * η_χ #* (F > -1e-8)
end

# special plastic helper functions
function compute_F(r::DruckerPrager, τ, P)
    F = (τ - P * r.sinϕ - r.C * r.cosϕ)
    return F*(F>-1e-8)
end
compute_Q(r::DruckerPrager, τ, P) = τ - P * r.sinψ

@inline function compute_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ - ε
    # return λ * ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) - ε
end

@inline function compute_volumetric_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return -λ * r.sinψ - θ
    # return λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) - θ
end

@inline compute_plastic_stress(r::DruckerPrager; τ_pl = 0, kwargs...) = τ_pl
@inline compute_stress(r::DruckerPrager; τ_pl = 0, kwargs...) = τ_pl

@inline compute_viscosity(r::DruckerPrager; kwargs...)   = Inf
@inline compute_viscosity_series(r::DruckerPrager; kwargs...)   = Inf
@inline compute_viscosity_parallel(r::DruckerPrager; kwargs...) = Inf
