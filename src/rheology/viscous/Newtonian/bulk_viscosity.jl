import ..RheologyCalculator: series_state_functions, parallel_state_functions, _isvolumetric
import ..RheologyCalculator: compute_volumetric_strain_rate, compute_pressure

"""
    BulkViscosity{T} <: AbstractViscosity

Represents the bulk viscosity of a material. Bulk viscosity is a material property that characterizes resistance to uniform compression or expansion.

# Fields
- `χ::T`: The value of the bulk viscosity.
"""
struct BulkViscosity{T} <: AbstractViscosity
    χ::T
end
@inline _isvolumetric(::BulkViscosity) = true
@inline series_state_functions(::BulkViscosity) = (compute_volumetric_strain_rate,)
@inline parallel_state_functions(::BulkViscosity) = (compute_pressure,)

@inline compute_volumetric_strain_rate(r::BulkViscosity; P = 0, kwargs...) = P / r.χ
@inline compute_pressure(r::BulkViscosity; θ = 0, kwargs...) = θ * r.χ
