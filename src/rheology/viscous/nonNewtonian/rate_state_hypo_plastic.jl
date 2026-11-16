# Hypoplastic Rate & State friction formulation
#
#module RateState_HypoPlastic

"""
    RateStateFriction{T, S} <: AbstractViscosity
    RateStateFriction(λ, μ₀, V₀, a, b, L, C, D; state_update = :explicit)

Implements rate and state dependent frictional viscous creeplas. 
Note that this is not a plasticity model, but a viscous model with state and rate dependence. 

`state_update` sets the state the friction law is evaluated at within a time step.
`:explicit` uses the state at the start of the step, `Ω_old`, and the state is
updated after the solve with `update_Ω`, as in Herrendörfer et al. (2018).
`:implicit` uses the state at the end of the step, `update_Ω(r; ε, Ω_old, dt)`,
inside the solve. The two schemes give different results.
`:implicit` is only available in a `ParallelModel`: in series the law gives the
strain rate from the stress, and with an implicit state that has no closed form.

# Fields
- `λ::T`    # fluid pressure ratio
- `μ₀::T`   # reference friction coefficient at reference slip velocity V₀
- `V₀::T`   # reference slip velocity
- `a::T`    # rate & state parameter    
- `b::T`    # rate & state parameter    
- `L::T`    # characteristic slip distance      
- `C::T`    # Cohesion
- `D::T`    # Cell size
"""
struct RateStateFriction{T, S} <: AbstractViscosity
    λ::T    # fluid pressure ratio
    μ₀::T   # reference friction coefficient at reference slip velocity V₀
    V₀::T   # reference slip velocity
    a::T    # rate & state parameter
    b::T    # rate & state parameter
    L::T    # characteristic slip distance
    C::T    # Cohesion
    D::T    # Cell size
end
function RateStateFriction(args...; state_update::Symbol = :explicit)
    state_update in (:explicit, :implicit) ||
        throw(ArgumentError("state_update must be :explicit or :implicit, got :$state_update"))
    p = promote(args...)
    return RateStateFriction{eltype(p), state_update}(p...)
end
@inline series_state_functions(::RateStateFriction) = (compute_strain_rate,)
check_series_leaf(::RateStateFriction{T, :implicit}) where {T} =
    throw(ArgumentError("RateStateFriction with state_update = :implicit only works in a ParallelModel"))
@inline parallel_state_functions(::RateStateFriction) = (compute_stress,)

"""
    state_in_step(r::RateStateFriction, ε, Ω_old, dt)

The state the friction law is evaluated at during a time step: `Ω_old` for
`state_update = :explicit`, the end-of-step state at strain rate `ε` for `:implicit`.
"""
@inline state_in_step(::RateStateFriction{T, :explicit}, ε, Ω_old, dt) where {T} = Ω_old
@inline state_in_step(r::RateStateFriction{T, :implicit}, ε, Ω_old, dt) where {T} = update_Ω(r; ε, Ω_old, dt)

# strain rate as a function of stress, state, and pressure
@inline function compute_strain_rate(r::RateStateFriction; τ = 0, Ω_old = 0, P = 0, dt = 0, kwargs...)

    Vp = 2 * r.V₀ * sinh(max((τ - r.C), 0) / (r.a * P * (1 - r.λ))) * exp(-(r.μ₀ + r.b * Ω_old) / r.a)

    ε = Vp / (2 * r.D)
    return ε
end

# stress as a function of strain rate, state, and pressure
@inline function compute_stress(r::RateStateFriction; ε = 0, Ω_old = 0, P = 0, dt = 0, kwargs...)
    Vp = 2 * r.D * ε
    Ω = state_in_step(r, ε, Ω_old, dt)
    μd = r.a * asinh(Vp / (2 * r.V₀) * exp((r.μ₀ + r.b * Ω) / r.a))
    τII = P * (1 - r.λ) * μd + r.C
    return τII
end

# This updates the state variable Ω based on slip velocity Vp and time step dt, which is needed outside the rheology definition.
# Θ = exp(Ω) is floored before the log, since it reaches zero for an overflowed slip rate or for dt = 0 with a fully relaxed Ω_old.
@inline function update_Ω(r::RateStateFriction; ε = 0, Ω_old = 0, dt = 0, kwargs...)
    Vp = 2 * r.D * ε
    if (Vp * dt / r.L ≤ 1.0e-6)
        Θ = exp(Ω_old) * (1 - Vp * dt / r.L) + r.V₀ * dt / r.L
    else
        Θ = r.V₀ / Vp + (exp(Ω_old) - r.V₀ / Vp) * exp(-Vp * dt / r.L)
    end
    return log(max(Θ, 1.0e-300))
end

# --------------------------------------------------------------------
