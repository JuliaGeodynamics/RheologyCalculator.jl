# Per-mechanism strain rates and dissipation.
#
# Deformation and heating are different questions and neither answers the other. Two
# elements can deform identically and dissipate by orders of magnitude differently
# (parallel branches, which share a strain rate but not a stress), or dissipate
# identically while deforming differently (series, which share a stress). Dissipation
# is therefore formed per element from that element's *own* conjugate pair, never by
# multiplying a global stress by a summed local rate.
#
# The pairing comes free from the equation structure. For a `compute_strain_rate`
# equation the per-leaf values are the elements' strain rates at the node's shared
# stress `args.τ`; for a `compute_stress` equation they are the elements' stresses at
# the node's shared rate `args.ε`. Same for the volumetric pair. Non-kinematic
# equations -- `compute_lambda` and friends, which are yield residuals rather than
# rates -- are skipped.
#
# Nothing here reads `others.τ0` or `others.P0`: no viscous or plastic state function
# takes history, and elastic elements are deliberately not reported, since elastic
# work is reversible storage and never enters a heat source.

"""
    ConstitutivePartition

Per-mechanism decomposition of the local deformation and its dissipation.

Rate fields are tuples and dissipation fields are scalars in W/m³, because the
energy equation wants one number per material point.

An entry is one *(kinematic equation, element)* contribution, not one element:
`viscous_ε[i]`, `viscous_θ[i]` and `viscous_mechanisms[i]` always refer to the same
contribution, with the mode an equation does not serve zeroed. An element appearing
in both a deviatoric and a volumetric equation -- as it does in any compressible
model -- therefore produces two entries, one carrying its `ε` and one its `θ`. The
scalars are unaffected: each mode is accumulated from its own equation only, so
nothing is double counted. Keeping all three tuples the same length is what lets a
single index mean the same thing in each; compacting to one entry per element would
require `ε` and `θ` to differ in length and break that.

The shear-heating source term is
`viscous_Φ + viscous_Φ_vol + plastic_Φ + plastic_Φ_vol`, optionally with a
Taylor-Quinney factor applied to the plastic pair. Elastic terms are absent by
design: elastic work is stored reversibly, not dissipated.
"""
struct ConstitutivePartition{Ev, Θv, Ep, Θp, Fv, Fvθ, Fp, Fpθ, Mv, Mp}
    viscous_ε::Ev          # deviatoric viscous strain rates
    viscous_θ::Θv          # volumetric viscous strain rates
    plastic_ε::Ep
    plastic_θ::Θp
    viscous_Φ::Fv          # Σ 2 τ_e ε_e
    viscous_Φ_vol::Fvθ     # Σ -P_e θ_e
    plastic_Φ::Fp
    plastic_Φ_vol::Fpθ
    viscous_mechanisms::Mv # tags aligned with viscous_ε / viscous_θ
    plastic_mechanisms::Mp
end

function Base.show(io::IO, ::MIME"text/plain", p::ConstitutivePartition)
    println(io, "ConstitutivePartition:")
    println(io, "  viscous  ε = ", p.viscous_ε, "  θ = ", p.viscous_θ)
    println(io, "           Φ = ", p.viscous_Φ, "  Φ_vol = ", p.viscous_Φ_vol)
    println(io, "  plastic  ε = ", p.plastic_ε, "  θ = ", p.plastic_θ)
    println(io, "           Φ = ", p.plastic_Φ, "  Φ_vol = ", p.plastic_Φ_vol)
    println(io, "  mechanisms: viscous ", p.viscous_mechanisms, ", plastic ", p.plastic_mechanisms)
    return nothing
end

"""
    shear_heating(p::ConstitutivePartition; β = 1)

Total irreversible dissipation, `W/m³`, suitable as a heat-equation source term.

`β` is the Taylor-Quinney factor applied to the plastic terms: experimentally a
fraction of plastic work is stored in microstructure rather than converted to heat,
typically `β ≈ 0.9`. There is no viscous counterpart for steady-state creep, whose
microstructure is stationary by definition; a dynamically recrystallising rheology
would need one.
"""
@inline shear_heating(p::ConstitutivePartition; β = 1) =
    p.viscous_Φ + p.viscous_Φ_vol + β * (p.plastic_Φ + p.plastic_Φ_vol)

# --- traversal ---------------------------------------------------------------

# Which conjugate an equation's per-leaf values pair with, and in which mode.
# `:dev_rate`  -> per-leaf are strain rates, conjugate is args.τ
# `:dev_stress`-> per-leaf are stresses,     conjugate is args.ε
# `:vol_rate`  -> per-leaf are volumetric rates, conjugate is args.P
# `:vol_press` -> per-leaf are pressures,    conjugate is args.θ
# anything else is not a conjugate pair and is skipped.
@inline _pair_mode(::F) where {F} = Val(:skip)
@inline _pair_mode(::typeof(compute_strain_rate)) = Val(:dev_rate)
@inline _pair_mode(::typeof(compute_stress)) = Val(:dev_stress)
@inline _pair_mode(::typeof(compute_volumetric_strain_rate)) = Val(:vol_rate)
@inline _pair_mode(::typeof(compute_pressure)) = Val(:vol_press)

@inline _mech_tag(r) = Val(nameof(typeof(r)))

# (deviatoric rate, volumetric rate, deviatoric Φ, volumetric Φ) for one leaf value
@inline _contribution(::Val{:dev_rate}, val, args) = (val, zero(val), 2 * args.τ * val, zero(val))
@inline _contribution(::Val{:dev_stress}, val, args) = (args.ε, zero(val), 2 * val * args.ε, zero(val))
@inline _contribution(::Val{:vol_rate}, val, args) = (zero(val), val, zero(val), -args.P * val)
@inline _contribution(::Val{:vol_press}, val, args) = (zero(val), args.θ, zero(val), -val * args.θ)

"""
    constitutive_partition(c::AbstractCompositeModel, sol::AbstractVector, vars, others)

Decompose the converged local solution into per-element strain rates and the
dissipation they produce. See [`ConstitutivePartition`](@ref).

`sol` must be a converged solution and `vars`/`others` the arguments that produced
it. `vars` is passed exactly as it was given to [`solve`](@ref); the elastic
correction of the prescribed strain rate is reapplied here so the partition is taken
at the state the Newton iteration actually solved.
"""
function constitutive_partition(c::AbstractCompositeModel, sol::AbstractVector, vars0, others)
    ε_corr = _direct_leaf_elastic_correction(c, vars0.ε, others)
    εII = second_invariant_value(vars0.ε .+ ε_corr)
    vars = merge(vars0, (; ε = εII))

    eqs = generate_equations(c)
    n = length(eqs)
    x = SVector{n}(ntuple(i -> @inbounds(sol[i]), n))
    args = generate_args_template(eqs, x, others)

    vε, vθ, vm = _gather(eqs, args, others, Val(:viscous))
    pε, pθ, pm = _gather(eqs, args, others, Val(:plastic))
    vΦ, vΦv = _gather_Φ(eqs, args, others, Val(:viscous))
    pΦ, pΦv = _gather_Φ(eqs, args, others, Val(:plastic))

    return ConstitutivePartition(vε, vθ, pε, pθ, vΦ, vΦv, pΦ, pΦv, vm, pm)
end

# Rates and tags for one category. Tuple recursion over equations and leaves; the
# category test is on the leaf type, so the selection folds at compile time.
# Fully unrolled, not recursive. The accumulator tuples grow as elements are
# visited, so a self-recursive call whose argument types widen at each level makes
# inference bail out to `Any` to guarantee termination. Unrolling gives every step
# its own SSA name and no recursion for inference to give up on.
@generated function _gather(eqs::NTuple{N, CompositeEquation}, args, others, ::Val{cat}) where {N, cat}
    N == 0 && return :(((), (), ()))
    return quote
        @inline
        ε_0 = (); θ_0 = (); m_0 = ()
        Base.@nexprs $N i -> begin
            (ε_i, θ_i, m_i) = _gather_eq(eqs[i], args[i], others, Val(cat), ε_{i - 1}, θ_{i - 1}, m_{i - 1})
        end
        ($(Symbol(:ε_, N)), $(Symbol(:θ_, N)), $(Symbol(:m_, N)))
    end
end

# Dispatch on the pair mode rather than branching on it: a non-kinematic equation
# returns the accumulator unchanged while a kinematic one returns extended tuples,
# so a runtime `if` would leave the return type unresolved.
@inline _gather_eq(eq::CompositeEquation, args, others, ::Val{cat}, ε, θ, m) where {cat} =
    _gather_eq(_pair_mode(eq.fn), eq, args, others, Val(cat), ε, θ, m)

@inline _gather_eq(::Val{:skip}, eq::CompositeEquation, args, others, ::Val{cat}, ε, θ, m) where {cat} =
    (ε, θ, m)

@inline function _gather_eq(mode::Val, eq::CompositeEquation, args, others, ::Val{cat}, ε, θ, m) where {cat}
    vals = evaluate_state_function_perleaf(eq.fn, eq.rheology, args, others, eq.el_number)
    return _accumulate_leaves(eq.rheology, vals, args, mode, Val(cat), ε, θ, m)
end

@generated function _accumulate_leaves(rheology::NTuple{NR, AbstractRheology}, vals, args, mode, ::Val{cat}, ε, θ, m) where {NR, cat}
    NR == 0 && return :((ε, θ, m))
    return quote
        @inline
        ε_0 = ε; θ_0 = θ; m_0 = m
        Base.@nexprs $NR i -> begin
            (ε_i, θ_i, m_i) = _leaf_step(rheology_category(rheology[i]), Val(cat), rheology[i], mode, vals[i], args, ε_{i - 1}, θ_{i - 1}, m_{i - 1})
        end
        ($(Symbol(:ε_, NR)), $(Symbol(:θ_, NR)), $(Symbol(:m_, NR)))
    end
end

@inline function _leaf_step(::Val{c}, ::Val{c}, r, mode, val, args, ε, θ, m) where {c}
    (dε, dθ, _, _) = _contribution(mode, val, args)
    return ((ε..., dε), (θ..., dθ), (m..., _mech_tag(r)))
end
@inline _leaf_step(::Val, ::Val, r, mode, val, args, ε, θ, m) = (ε, θ, m)


# Dissipation scalars for one category.
@generated function _gather_Φ(eqs::NTuple{N, CompositeEquation}, args, others, ::Val{cat}) where {N, cat}
    N == 0 && return :((0.0, 0.0))
    return quote
        @inline
        Φ_0 = 0.0; Φv_0 = 0.0
        Base.@nexprs $N i -> begin
            (Φ_i, Φv_i) = _gather_Φ_eq(eqs[i], args[i], others, Val(cat), Φ_{i - 1}, Φv_{i - 1})
        end
        ($(Symbol(:Φ_, N)), $(Symbol(:Φv_, N)))
    end
end

@inline _gather_Φ_eq(eq::CompositeEquation, args, others, ::Val{cat}, Φ, Φv) where {cat} =
    _gather_Φ_eq(_pair_mode(eq.fn), eq, args, others, Val(cat), Φ, Φv)

@inline _gather_Φ_eq(::Val{:skip}, eq::CompositeEquation, args, others, ::Val{cat}, Φ, Φv) where {cat} =
    (Φ, Φv)

@inline function _gather_Φ_eq(mode::Val, eq::CompositeEquation, args, others, ::Val{cat}, Φ, Φv) where {cat}
    vals = evaluate_state_function_perleaf(eq.fn, eq.rheology, args, others, eq.el_number)
    return _accumulate_Φ(eq.rheology, vals, args, mode, Val(cat), Φ, Φv)
end

@generated function _accumulate_Φ(rheology::NTuple{NR, AbstractRheology}, vals, args, mode, ::Val{cat}, Φ, Φv) where {NR, cat}
    NR == 0 && return :((Φ, Φv))
    return quote
        @inline
        Φ_0 = Φ; Φv_0 = Φv
        Base.@nexprs $NR i -> begin
            (Φ_i, Φv_i) = _Φ_step(rheology_category(rheology[i]), Val(cat), mode, vals[i], args, Φ_{i - 1}, Φv_{i - 1})
        end
        ($(Symbol(:Φ_, NR)), $(Symbol(:Φv_, NR)))
    end
end

@inline function _Φ_step(::Val{c}, ::Val{c}, mode, val, args, Φ, Φv) where {c}
    (_, _, dΦ, dΦv) = _contribution(mode, val, args)
    return (Φ + dΦ, Φv + dΦv)
end
@inline _Φ_step(::Val, ::Val, mode, val, args, Φ, Φv) = (Φ, Φv)
