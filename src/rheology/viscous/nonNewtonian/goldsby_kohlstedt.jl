"""Goldsby–Kohlstedt ice creep mechanisms.

The full ice law is assembled with `ParallelModel` and `SeriesModel`; these
types provide only the mechanisms not covered by the generic creep laws.
"""

struct GoldsbyKohlstedtDiffusion{T} <: AbstractViscosity
    R::T
    D0v::T
    Qv::T
    D0b::T
    Qb::T
    Vm::T
    δ::T
end
GoldsbyKohlstedtDiffusion(args...) = GoldsbyKohlstedtDiffusion(promote(args...)...)

@inline series_state_functions(::GoldsbyKohlstedtDiffusion) = (compute_strain_rate,)
@inline parallel_state_functions(::GoldsbyKohlstedtDiffusion) = (compute_stress,)

@inline function _gk_diffusion_viscosity(r, T, d)
    (; R, D0v, Qv, D0b, Qb, Vm, δ) = r
    D = D0v * exp(-Qv / (R * T)) + π * δ / d * D0b * exp(-Qb / (R * T))
    return R * T * d^2 / (126 * Vm * D)
end

@inline compute_strain_rate(r::GoldsbyKohlstedtDiffusion; τ = 0, T = 0, d = 1, kwargs...) =
    τ / (2 * _gk_diffusion_viscosity(r, T, d))
@inline compute_stress(r::GoldsbyKohlstedtDiffusion; ε = 0, T = 0, d = 1, kwargs...) =
    2 * _gk_diffusion_viscosity(r, T, d) * ε
@inline compute_viscosity(r::GoldsbyKohlstedtDiffusion; T = 0, d = 1, kwargs...) = _gk_diffusion_viscosity(r, T, d)

struct GoldsbyKohlstedtCreep{N, P, T} <: AbstractViscosity
    n::N
    p::P
    A::T
    Q::T
    R::T
end
GoldsbyKohlstedtCreep(n, p, args...) = GoldsbyKohlstedtCreep(n, p, promote(args...)...)

@inline series_state_functions(::GoldsbyKohlstedtCreep) = (compute_strain_rate,)
@inline parallel_state_functions(::GoldsbyKohlstedtCreep) = (compute_stress,)

@inline function compute_strain_rate(r::GoldsbyKohlstedtCreep; τ = 0, T = 0, d = 1, kwargs...)
    (; n, p, A, Q, R) = r
    return A * 3^((n + 1) / 2) * τ^n * d^(-p) * exp(-Q / (R * T)) / 2
end

@inline function compute_stress(r::GoldsbyKohlstedtCreep; ε = 0, T = 0, d = 1, kwargs...)
    (; n, p, A, Q, R) = r
    return (2 * abs(ε) * d^p / (A * 3^((n + 1) / 2)) * exp(Q / (R * T)))^(1 / n) * sign(ε)
end
