"""
    IncompressibleElasticity{T} <: AbstractElasticity

Represents incompressible elastic deformation (shear only).

# Fields
- `G::T`: The shear modulus.
"""
struct IncompressibleElasticity{T} <: AbstractElasticity
    G::T
end
@inline series_state_functions(::IncompressibleElasticity) = (compute_strain_rate,)
@inline parallel_state_functions(::IncompressibleElasticity) = (compute_stress,)

@inline compute_strain_rate(r::IncompressibleElasticity; τ = 0, τ0 = 0, dt = 0, kwargs...) = τ / (2 * r.G * dt)
@inline compute_stress(r::IncompressibleElasticity; ε = 0, τ0 = 0, dt = 0, kwargs...) = 2 * r.G * dt * ε

@inline compute_viscosity(r::IncompressibleElasticity; dt = 0, kwargs...) = r.G * dt
