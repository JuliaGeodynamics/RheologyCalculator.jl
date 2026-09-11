# Ice rheology following Goldsby & Kohlstedt (2001) as formulated in
# Kihoulou et al. (2025), Sci. Adv. 11, eadq8719.
#
# The full ductile ice viscosity (Eq. 5) combines four mechanisms in parallel:
#
#   1/η_duct = 1/η_diff + 1/η_disl + 1/(η_GBS + η_BS) + 1/η_max
#
# where GBS and BS are in series (their viscosities add). Using this package:
#
#   diff    = GoldsbyKohlstedtDiffusion(...)
#   disl    = GoldsbyKohlstedtCreep(n_disl, 0, A_disl, Q_disl, R)
#   gbs     = GoldsbyKohlstedtCreep(n_gbs,  p_gbs,  A_gbs,  Q_gbs,  R)
#   bs      = GoldsbyKohlstedtCreep(n_bs,   p_bs,   A_bs,   Q_bs,   R)
#   η_max   = LinearViscosity(η_max)
#
#   c = SeriesModel(ParallelModel(diff, disl, SeriesModel(gbs, bs), η_max))
#
# Material parameters for H₂O ice Ih are listed in Table S2 of Kihoulou et al.
# (originally from Goldsby & Kohlstedt 2001, J. Geophys. Res.).

using RheologyCalculator
import RheologyCalculator: series_state_functions, parallel_state_functions
import RheologyCalculator: compute_strain_rate, compute_stress
import RheologyCalculator: compute_viscosity, compute_viscosity_series, compute_viscosity_parallel

# ── GoldsbyKohlstedtDiffusion ──────────────────────────────────────────────

"""
    GoldsbyKohlstedtDiffusion{T} <: AbstractViscosity

Ice diffusion creep combining volume and grain-boundary diffusion channels
(Goldsby & Kohlstedt 2001, Eq. 6 of Kihoulou et al. 2025):

    η_diff = R T d² / (3/2 · 84 Vₘ · [D₀ᵥ exp(-Qᵥ/RT) + π δ/d · D₀ᵦ exp(-Qᵦ/RT)])

# Fields
- `R::T`: universal gas constant (J mol⁻¹ K⁻¹)
- `D0v::T`: volume-diffusion pre-exponential (m² s⁻¹)
- `Qv::T`: volume-diffusion activation energy (J mol⁻¹)
- `D0b::T`: grain-boundary-diffusion pre-exponential (m³ s⁻¹)
- `Qb::T`: grain-boundary-diffusion activation energy (J mol⁻¹)
- `Vm::T`: molar volume of ice (m³ mol⁻¹)
- `δ::T`: grain-boundary width (m)

# Keyword arguments consumed at call time (from `others`)
- `T`: temperature (K)
- `d`: grain size (m) — treated as a per-element history field
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

@inline function _gk_diff_viscosity(r::GoldsbyKohlstedtDiffusion, T, d)
    (; R, D0v, Qv, D0b, Qb, Vm, δ) = r
    Deff = D0v * exp(-Qv / (R * T)) + π * δ / d * D0b * exp(-Qb / (R * T))
    return R * T * d^2 / (126 * Vm * Deff)  # 3/2 * 84 = 126
end

@inline compute_strain_rate(r::GoldsbyKohlstedtDiffusion; τ = 0, T = 0, d = 1, kwargs...) =
    τ / (2 * _gk_diff_viscosity(r, T, d))

@inline compute_stress(r::GoldsbyKohlstedtDiffusion; ε = 0, T = 0, d = 1, kwargs...) =
    2 * _gk_diff_viscosity(r, T, d) * ε

@inline compute_viscosity(r::GoldsbyKohlstedtDiffusion; T = 0, d = 1, kwargs...) = _gk_diff_viscosity(r, T, d)
@inline compute_viscosity_series(r::GoldsbyKohlstedtDiffusion; T = 0, d = 1, kwargs...) = _gk_diff_viscosity(r, T, d)
@inline compute_viscosity_parallel(r::GoldsbyKohlstedtDiffusion; T = 0, d = 1, kwargs...) = _gk_diff_viscosity(r, T, d)

# ── GoldsbyKohlstedtCreep ──────────────────────────────────────────────────

"""
    GoldsbyKohlstedtCreep{I, T} <: AbstractViscosity

Stress-dependent ice creep (dislocation creep, basal slip, or grain-boundary
sliding) following Goldsby & Kohlstedt (2001), Eq. 7 of Kihoulou et al. (2025).

The effective viscosity is:

    η = d^p / (3^((n+1)/2) · A · σII^(n-1)) · exp(Q / RT)

giving the strain rate:

    ε = A · 3^((n+1)/2) · τ^n · d^(-p) / 2 · exp(-Q / RT)

Instantiate once per mechanism with its own (n, p, A, Q) and compose with
`ParallelModel` / `SeriesModel` as needed (see module docstring above).

# Fields
- `n::I`: stress exponent
- `p::I`: grain-size exponent (0 for dislocation creep)
- `A::T`: pre-exponential factor (Pa⁻ⁿ mᵖ s⁻¹)
- `Q::T`: activation energy (J mol⁻¹)
- `R::T`: universal gas constant (J mol⁻¹ K⁻¹)

# Keyword arguments consumed at call time (from `others`)
- `T`: temperature (K)
- `d`: grain size (m) — treated as a per-element history field
"""
struct GoldsbyKohlstedtCreep{I, T} <: AbstractViscosity
    n::I
    p::I
    A::T
    Q::T
    R::T
end
GoldsbyKohlstedtCreep(n, p, args...) = GoldsbyKohlstedtCreep(n, p, promote(args...)...)

@inline series_state_functions(::GoldsbyKohlstedtCreep) = (compute_strain_rate,)
@inline parallel_state_functions(::GoldsbyKohlstedtCreep) = (compute_stress,)

@inline function _gk_creep_viscosity(r::GoldsbyKohlstedtCreep, τ, T, d)
    (; n, p, A, Q, R) = r
    τ_eff = max(abs(τ), eps(typeof(float(τ))))
    return d^p / (A * 3^((n + 1) / 2) * τ_eff^(n - 1)) * exp(Q / (R * T))
end

@inline function compute_strain_rate(r::GoldsbyKohlstedtCreep; τ = 0, T = 0, d = 1, kwargs...)
    η = _gk_creep_viscosity(r, τ, T, d)
    return τ / (2 * η)
end

@inline function compute_stress(r::GoldsbyKohlstedtCreep; ε = 0, T = 0, d = 1, kwargs...)
    (; n, p, A, Q, R) = r
    return (2 * abs(ε) * d^p / (A * 3^((n + 1) / 2)) * exp(Q / (R * T)))^(1 / n) * sign(ε)
end

@inline compute_viscosity(r::GoldsbyKohlstedtCreep; τ = 0, T = 0, d = 1, kwargs...) =
    _gk_creep_viscosity(r, τ, T, d)

@inline compute_viscosity_series(r::GoldsbyKohlstedtCreep; ε = 0, T = 0, d = 1, kwargs...) =
    compute_stress(r; ε = ε, T = T, d = d) / (2 * max(abs(ε), eps(typeof(float(ε)))))

@inline compute_viscosity_parallel(r::GoldsbyKohlstedtCreep; τ = 0, T = 0, d = 1, kwargs...) =
    _gk_creep_viscosity(r, τ, T, d)
