"""
    LinearViscosity{T} <: AbstractViscosity

Represents a linear viscosity model following Newton's law of viscosity.

# Fields
- `η::T`: The dynamic viscosity coefficient.
"""
struct LinearViscosity{T} <: AbstractViscosity
    η::T
end
@inline series_state_functions(::LinearViscosity) = (compute_strain_rate,)
@inline parallel_state_functions(::LinearViscosity) = (compute_stress,)

@inline compute_strain_rate(r::LinearViscosity; τ = 0, kwargs...) = τ / (2 * r.η)
@inline compute_stress(r::LinearViscosity; ε = 0, kwargs...) = ε * 2 * r.η
@inline compute_viscosity(r::LinearViscosity; kwargs...) = r.η
