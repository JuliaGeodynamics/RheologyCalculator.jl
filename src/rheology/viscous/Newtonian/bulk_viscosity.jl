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

# Sign convention: `P` is positive in compression and `θ` positive in dilation, so
# compressing the material must compact it and the two carry opposite signs. This
# matches `BulkElasticity`, whose `θ = -(P - P0)/(K*dt)` carries the same minus.
# Without it the volumetric dissipation `-P*θ = -P^2/χ` would be negative, i.e. a
# viscous element generating energy.
@inline compute_volumetric_strain_rate(r::BulkViscosity; P = 0, kwargs...) = -P / r.χ
@inline compute_pressure(r::BulkViscosity; θ = 0, kwargs...) = -θ * r.χ
