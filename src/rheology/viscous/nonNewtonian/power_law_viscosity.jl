import ..RheologyCalculator: series_state_functions, parallel_state_functions
import ..RheologyCalculator: compute_strain_rate, compute_stress
import ..RheologyCalculator: compute_viscosity, compute_viscosity_series, compute_viscosity_parallel

"""
    PowerLawViscosity{T,I} <: AbstractViscosity

Represents a power-law viscosity model where viscosity depends on strain rate.

# Fields
- `η::T`: The viscosity coefficient.
- `n::I`: The power-law exponent (note: not promoted to floating point by default).
"""
struct PowerLawViscosity{T, I} <: AbstractViscosity
    η::T
    n::I # DO NOT PROMOTE TO FP BY DEFAULT
end
@inline series_state_functions(::PowerLawViscosity) = (compute_strain_rate,)
@inline parallel_state_functions(::PowerLawViscosity) = (compute_stress,)

@inline compute_strain_rate(r::PowerLawViscosity; τ = 0, kwargs...) = τ^r.n / (2 * r.η)
@inline compute_stress(r::PowerLawViscosity; ε = 0, kwargs...) = ε^(1 / r.n) * (2 * r.η)^(1 / r.n)

@inline compute_viscosity(r::PowerLawViscosity; ε = 0, kwargs...)   = compute_stress(r; ε = ε, kwargs...)/(2*ε)
@inline compute_viscosity_series(r::PowerLawViscosity; ε = 0, kwargs...)   = compute_stress(r; ε = ε, kwargs...)/(2*ε)
@inline compute_viscosity_parallel(r::PowerLawViscosity; τ = 0, kwargs...) = τ / (2 * compute_strain_rate(r; τ = τ, kwargs...))
