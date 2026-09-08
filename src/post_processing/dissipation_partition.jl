
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

function Base.show(io::IO, ::MIME"text/plain", p::DissipationPartition)
    println(io, "DissipationPartition:")
    println(io, "  viscous  Φ = ", p.viscous_Φ)
    println(io, "  plastic  Φ = ", p.plastic_Φ)
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
`vars` is accepted to preserve the conventional package calling sequence; unlike
deformation-rate diagnostics, dissipation is reconstructed entirely from `sol` and
`others`.
"""
function dissipation_partition(c::AbstractCompositeModel, sol::AbstractVector, _vars, others)
    eqs = generate_equations(c)
    n = length(eqs)
    x = SVector{n}(ntuple(i -> @inbounds(sol[i]), n))
    args = generate_args_template(eqs, x, others)

    vΦ, pΦ = _dissipation(eqs, args, others)
    return DissipationPartition(vΦ, pΦ)
end

# Walk all equations once. Generated unrolling keeps the traversal type-stable and
# allocation-free for the statically shaped equation tuple.
@generated function _dissipation(eqs::NTuple{N, CompositeEquation}, args, others) where {N}
    N == 0 && return :((0.0, 0.0))
    return quote
        @inline
        viscous_0 = 0.0
        plastic_0 = 0.0
        Base.@nexprs $N i -> begin
            viscous_i, plastic_i = _equation_dissipation(
                eqs[i].fn,
                eqs[i],
                args[i],
                others,
                viscous_{i - 1},
                plastic_{i - 1},
            )
        end
        ($(Symbol(:viscous_, N)), $(Symbol(:plastic_, N)))
    end
end

# Only deviatoric kinematic equations form stress-strain-rate conjugate pairs.
# Yield, volumetric, and other equations do not contribute to shear heating
# and are therefore ignored to be more efficient.
@inline _equation_dissipation(::F, _eq, _args, _others, viscous, plastic) where {F} =
    (viscous, plastic)

# here, we evaluate because we know that there is viscous or plastic deviatoric
@inline function _equation_dissipation(
        fn::Union{typeof(compute_strain_rate), typeof(compute_stress)},
        eq,
        args,
        others,
        viscous,
        plastic,
    )
    values = evaluate_state_function_perleaf(fn, eq.rheology, args, others, eq.el_number)
    return _leaf_dissipation(eq.rheology, values, fn, args, viscous, plastic)
end

@generated function _leaf_dissipation(
        rheology::NTuple{N, AbstractRheology},
        values,
        fn,
        args,
        viscous,
        plastic,
    ) where {N}
    N == 0 && return :((viscous, plastic))
    return quote
        @inline
        viscous_0 = viscous
        plastic_0 = plastic
        Base.@nexprs $N i -> begin
            viscous_i, plastic_i = _add_dissipation(
                rheology_category(rheology[i]),
                fn,
                values[i],
                args,
                viscous_{i - 1},
                plastic_{i - 1},
            )
        end
        ($(Symbol(:viscous_, N)), $(Symbol(:plastic_, N)))
    end
end

@inline _deviatoric_power(::typeof(compute_strain_rate), value, args) =
    2 * args.τ * value
@inline _deviatoric_power(::typeof(compute_stress), value, args) =
    2 * value * args.ε

@inline _add_dissipation(::Val{:viscous}, fn, value, args, viscous, plastic) =
    (viscous + _deviatoric_power(fn, value, args), plastic)
@inline _add_dissipation(::Val{:plastic}, fn, value, args, viscous, plastic) =
    (viscous, plastic + _deviatoric_power(fn, value, args))
@inline _add_dissipation(::Val, fn, value, args, viscous, plastic) =
    (viscous, plastic)
