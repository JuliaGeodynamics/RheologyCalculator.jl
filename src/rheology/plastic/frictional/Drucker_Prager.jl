"""
    DruckerPrager(C, ϕ, ψ)
    DruckerPrager(C, ϕ, ψ, η_vp)

Represents the Drucker-Prager plasticity model for pressure-dependent yielding.
ϕ and ψ are the friction and dilatancy angles in degrees. The optional
Duvaut-Lions regularisation viscosity `η_vp` uses the same viscosity units as
the rest of the model. It defaults to `one(C)` for backwards compatibility;
set it to zero to recover the unregularised yield condition.

The second type parameter records whether the element is dilatant, so that the
volumetric equations a dilatant element needs are added to the local system at
compile time and cost nothing when `ψ` is zero.

# Fields
- `C::T`: The cohesion parameter.
- `ϕ::T`: The friction angle (in degrees).
- `ψ::T`: The dilatancy angle (in degrees).
- `η_vp::T`: The Duvaut-Lions viscoplastic regularisation viscosity.
"""
struct DruckerPrager{T, D} <: AbstractPlasticity
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
        return new{T, !iszero(sinψ)}(C, ϕ, ψ, η_vp, sinϕ, sinψ, cosϕ, cosψ)
    end
end
function DruckerPrager(C, ϕ, ψ)
    C, ϕ, ψ = promote(C, ϕ, ψ)
    return DruckerPrager(C, ϕ, ψ, one(C))
end
DruckerPrager(args::Vararg{Any, 4}) = DruckerPrager(promote(args...)...)

@inline _isvolumetric(::DruckerPrager{T, D}) where {T, D} = D

@inline series_state_functions(::DruckerPrager{T, false}) where {T} = (compute_strain_rate, compute_lambda)
@inline series_state_functions(::DruckerPrager{T, true}) where {T} = (compute_strain_rate, compute_volumetric_strain_rate, compute_lambda)

@inline parallel_state_functions(::DruckerPrager{T, false}) where {T} = (compute_stress, compute_pressure, compute_lambda_parallel, compute_plastic_strain_rate)
@inline parallel_state_functions(::DruckerPrager{T, true}) where {T} = (compute_stress, compute_pressure, compute_lambda_parallel, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate)

@inline function compute_strain_rate(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r::DruckerPrager; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    return ε_pl / 2
end

# Dilatant flow rule, θ_pl = -λ ∂Q/∂P = λ sinψ. `θ` is positive in dilation, so a
# dilatant element expands as it shears. The deviatoric plastic strain-rate
# invariant is λ/2, hence the factor of two between them: θ_pl = 2 sinψ ε̇_II_pl.
@inline compute_volumetric_strain_rate(r::DruckerPrager; λ = 0, kwargs...) = λ * r.sinψ

@inline function compute_lambda(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    return F - λ * r.η_vp
end

@inline function compute_lambda_parallel(r::DruckerPrager; τ_pl = 0, λ = 0, P = 0, P_pl = 0, kwargs...)
    P_yield = _parallel_pressure(r, P, P_pl)
    F = compute_F(r, τ_pl, P_yield)
    return F - λ * r.η_vp
end

# The pressure the yield surface sees inside a parallel composite. A dilatant
# element carries its own pressure unknown, shared with the rest of the branch;
# a non-dilatant one has none, and `P` is then the composite's pressure.
@inline _parallel_pressure(::DruckerPrager{T, false}, P, P_pl) where {T} = P
@inline _parallel_pressure(::DruckerPrager{T, true}, P, P_pl) where {T} = P_pl

# special plastic helper functions
function compute_F(r::DruckerPrager, τ, P)
    F = (τ - P * r.sinϕ - r.C * r.cosϕ)
    return F * (F > -1.0e-8)
end
compute_Q(r::DruckerPrager, τ, P) = τ - P * r.sinψ

@inline function compute_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ - ε
end

# The parallel counterpart of `compute_volumetric_strain_rate`: every element of a
# parallel composite shares the branch's volumetric strain rate, which the residual
# assembly subtracts, so the two carry the same expression.
@inline compute_volumetric_plastic_strain_rate(r::DruckerPrager; λ = 0, kwargs...) = λ * r.sinψ

@inline compute_plastic_stress(r::DruckerPrager; τ_pl = 0, kwargs...) = τ_pl
@inline compute_stress(r::DruckerPrager; τ_pl = 0, kwargs...) = τ_pl
@inline compute_pressure(r::DruckerPrager; P_pl = 0, kwargs...) = P_pl

@inline compute_viscosity(r::DruckerPrager; kwargs...) = Inf
