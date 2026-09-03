import ..RheologyCalculator: series_state_functions, parallel_state_functions, _isvolumetric
import ..RheologyCalculator: compute_strain_rate, compute_stress, compute_pressure, compute_volumetric_strain_rate
import ..RheologyCalculator: compute_viscosity, compute_viscosity_series, compute_viscosity_parallel

"""
    Elasticity{T} <: AbstractElasticity

Represents elastic deformation with both shear and bulk components.

# Fields
- `G::T`: The shear modulus.
- `K::T`: The bulk modulus.
"""
struct Elasticity{T} <: AbstractElasticity
    G::T
    K::T
end
@inline _isvolumetric(::Elasticity) = true
@inline series_state_functions(::Elasticity) = (compute_strain_rate, compute_volumetric_strain_rate)
@inline parallel_state_functions(::Elasticity) = (compute_stress, compute_pressure)

@inline compute_strain_rate(r::Elasticity; τ = 0, τ0 = 0, dt = 0, kwargs...) = τ / (2 * r.G * dt)
@inline compute_volumetric_strain_rate(r::Elasticity; P = 0, P0 = 0, dt = 0, kwargs...) = -(P - P0) / (r.K * dt)
@inline compute_stress(r::Elasticity; ε = 0, τ0 = 0, dt = 0, kwargs...) = 2 * r.G * dt * ε
@inline compute_pressure(r::Elasticity; θ = 0, P0 = 0, dt = 0, kwargs...) = P0 - r.K * dt * θ

@inline compute_viscosity(r::Elasticity; dt = 0, kwargs...)   = r.G * dt
@inline compute_viscosity_series(r::Elasticity; dt = 0, kwargs...)   = r.G * dt
@inline compute_viscosity_parallel(r::Elasticity; dt = 0, kwargs...) = r.G * dt
