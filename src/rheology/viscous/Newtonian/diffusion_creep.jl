import ..RheologyCalculator: series_state_functions, parallel_state_functions
import ..RheologyCalculator: compute_strain_rate, compute_stress
import ..RheologyCalculator: compute_viscosity, compute_viscosity_series, compute_viscosity_parallel

"""
    DiffusionCreep{I,T} <: AbstractViscosity

Represents diffusion creep deformation mechanism in materials.

# Fields
- `n::I`: The stress exponent.
- `r::T`: The water fugacity exponent.
- `p::T`: The grain size exponent.
- `A::T`: The material-specific rheological parameter.
- `E::T`: The activation energy.
- `V::T`: The activation volume.
- `R::T`: The universal gas constant.
"""
struct DiffusionCreep{I, T} <: AbstractViscosity
    n::I
    r::T
    p::T
    A::T
    E::T
    V::T
    R::T
end
DiffusionCreep(args...) = DiffusionCreep(args[1], promote(args[2:end]...)...)
@inline series_state_functions(::DiffusionCreep) = (compute_strain_rate,)
@inline parallel_state_functions(::DiffusionCreep) = (compute_stress,)

@inline function compute_strain_rate(r::DiffusionCreep; τ = 0, T = 0, P = 0, f = 1, d = 1, args...)
    (; n, r, p, A, E, V, R) = r

    ε = A * τ^n * f^r * d^(-p) * exp(-(E + P * V) / (R * T))
    return ε
end

@inline function compute_stress(r::DiffusionCreep; ε = 0, T = 0, P = 0, f = 1, d = 1, args...)
    (; n, r, p, A, E, V, R) = r

    _n = inv(n)
    τ = A^(-_n) * ε^_n * f^(-r * _n) * d^(p * _n) * exp((E + P * V) / (n * R * T))

    return τ
end

@inline compute_viscosity(r::DiffusionCreep; ε = 0, kwargs...)   = compute_stress(r; ε = ε, kwargs...)/(2*ε)
@inline compute_viscosity_series(r::DiffusionCreep; ε = 0, kwargs...)   = compute_stress(r; ε = ε, kwargs...)/(2*ε)
@inline compute_viscosity_parallel(r::DiffusionCreep; τ = 0, kwargs...) = τ / (2 * compute_strain_rate(r; τ = τ, kwargs...))
