import ..RheologyCalculator: series_state_functions
import ..RheologyCalculator: compute_strain_rate, compute_stress
import ..RheologyCalculator: compute_viscosity, compute_viscosity_series, compute_viscosity_parallel

"""
    LTPViscosity{T} <: AbstractViscosity

Represents a low-temperature plasticity (LTP) viscosity model.

# Fields
- `ε0::T`: The reference strain rate (default: 6.2e-13).
- `Q::T`: The activation energy (default: 76).
- `σb::T`: The brittle strength (default: 1.8 GPa).
- `σr::T`: The reference stress (default: 3.4 GPa).
"""
struct LTPViscosity{T} <: AbstractViscosity
    ε0::T # 6.2e-13
    Q::T  # 76
    σb::T # 1.8 GPa
    σr::T # 3.4 Gpa
end
LTPViscosity(args...) = LTPViscosity(promote(args...)...)
@inline series_state_functions(::LTPViscosity) = (compute_strain_rate,)

@inline compute_strain_rate(r::LTPViscosity; τ = 0e0, kwargs...) = max(r.ε0 * sinh(r.Q * (τ - r.σb) / r.σr), 0.0)
# @inline compute_strain_rate(r::LTPViscosity; τ = 0, kwargs...) = r.ε0 * sinh(r.Q * (τ - r.σb) / r.σr)
@inline compute_stress(r::LTPViscosity; ε = 0e0, kwargs...) = r.σr / r.Q * asinh(ε / r.ε0) + r.σb

@inline compute_viscosity(r::LTPViscosity; ε = 0e0, kwargs...)   = compute_stress(r; ε = ε, kwargs...)/(2*ε)
@inline compute_viscosity_series(r::LTPViscosity; ε = 0e0, kwargs...)   = compute_stress(r; ε = ε, kwargs...)/(2*ε)
@inline compute_viscosity_parallel(r::LTPViscosity; τ = 00e, kwargs...) = τ / (2 * compute_strain_rate(r; τ = τ, kwargs...))
