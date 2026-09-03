import ..RheologyCalculator: series_state_functions, parallel_state_functions, _isvolumetric
import ..RheologyCalculator: compute_volumetric_strain_rate, compute_pressure

"""
    BulkElasticity{T} <: AbstractElasticity

Represents bulk elastic deformation (volumetric compression/expansion only).

# Fields
- `K::T`: The bulk modulus.
"""
struct BulkElasticity{T} <: AbstractElasticity
    K::T
end
@inline _isvolumetric(::BulkElasticity) = true
@inline series_state_functions(::BulkElasticity) = (compute_volumetric_strain_rate,)
@inline parallel_state_functions(::BulkElasticity) = (compute_pressure,)
@inline compute_volumetric_strain_rate(r::BulkElasticity; P = 0, P0 = 0, dt = 0, kwargs...) = -(P - P0) / (r.K * dt)
@inline compute_pressure(r::BulkElasticity; θ = 0, P0 = 0, dt = 0, kwargs...) = P0 - r.K * dt * θ
