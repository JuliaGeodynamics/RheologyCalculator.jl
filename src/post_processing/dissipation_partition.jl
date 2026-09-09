
"""
    DissipationPartition

Decomposition of irreversible deviatoric dissipation by rheological category.

Both fields are scalars in W/m³, ready to be used in an energy equation. Elastic
terms are absent by design: elastic work is reversible storage and is not
dissipated.
"""
struct DissipationPartition{Fv, Fp}
    viscous_Φ::Fv  # total deviatoric viscous dissipation
    plastic_Φ::Fp  # total deviatoric plastic dissipation
end

# `_show_round` is defined in component_partition.jl, which is included first.

function Base.show(io::IO, ::MIME"text/plain", p::DissipationPartition)
    println(io, "DissipationPartition:")
    println(io, "  viscous  Φ = ", _show_round(p.viscous_Φ))
    println(io, "  plastic  Φ = ", _show_round(p.plastic_Φ))
    return nothing
end

"""
    shear_heating(p::DissipationPartition; β = 1)

Deviatoric irreversible dissipation, `W/m³`, suitable as a shear-heating
source term.

`β` is the Taylor-Quinney factor applied to the plastic terms: experimentally a
fraction of plastic work is stored in microstructure rather than converted to heat,
typically `β ≈ 0.9`. There is no viscous counterpart for steady-state creep, whose
microstructure is stationary.
"""
@inline shear_heating(p::DissipationPartition; β = 1) =
    p.viscous_Φ + β * p.plastic_Φ

"""
    dissipation_partition(c::AbstractCompositeModel, sol::AbstractVector, vars, others)

Decompose the converged local solution into viscous and plastic dissipation. See
[`DissipationPartition`](@ref).

`sol`, `vars`, and `others` must be the arguments supplied to [`solve`](@ref).
"""
function dissipation_partition(c::AbstractCompositeModel, sol::AbstractVector, vars, others)
    p = component_partition(c, sol, vars, others)
    return DissipationPartition(_total_power(p.viscous_τ, p.viscous_ε), _total_power(p.plastic_τ, p.plastic_ε))
end

@inline _total_power(τ::Tuple{}, ε::Tuple{}) = 0.0
@inline _total_power(τ::Tuple, ε::Tuple) =
    sum(2 .* τ .* ε; init = zero(eltype(τ)) * zero(eltype(ε)))
