import ..RheologyCalculator: series_state_functions, parallel_state_functions
import ..RheologyCalculator: compute_strain_rate, compute_stress

"""
    DislocationCreep{I,T} <: AbstractViscosity

Represents dislocation creep deformation mechanism in materials.

# Fields
- `n::I`: The power-law exponent.
- `r::T`: The exponent of water-fugacity.
- `A::T`: The material specific rheological parameter.
- `E::T`: The activation energy.
- `V::T`: The activation volume.
- `R::T`: The universal gas constant.
"""
struct DislocationCreep{I, T} <: AbstractViscosity
    n::I # power-law exponent
    r::T # exponent of water-fugacity
    A::T # material specific rheological parameter
    E::T # activation energy
    V::T # activation volume
    R::T # universal gas constant
end
DislocationCreep(args...) = DislocationCreep(args[1], promote(args[2:end]...)...)
@inline series_state_functions(::DislocationCreep) = (compute_strain_rate,)
@inline parallel_state_functions(::DislocationCreep) = (compute_stress,)

@inline function compute_strain_rate(r::DislocationCreep; τ = 0, T = 0, P = 0, f = 0, args...)
    (; n, r, A, E, V, R) = r

    ε = A * τ^n * f^r * exp(-(E + P * V) / (R * T))
    return ε
end

@inline function compute_stress(r::DislocationCreep; ε = 0, T = 0, P = 0, f = 0, args...)
    (; n, r, A, E, V, R) = r

    _n = inv(n)

    τ = A^(-_n) * ε^_n * f^(-r * _n) * exp((E + P * V) / (n * R * T))

    return τ
end
