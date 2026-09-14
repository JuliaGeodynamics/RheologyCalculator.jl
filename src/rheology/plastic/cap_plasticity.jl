# State functions shared by every AbstractCapPlasticity element. A subtype
# supplies only its struct, compute_F, and compute_Q; everything below is
# expressed in terms of those two and so needs no per-model definition.

@inline _isvolumetric(::AbstractCapPlasticity) = true

@inline series_state_functions(::AbstractCapPlasticity) = (compute_strain_rate, compute_lambda, compute_volumetric_strain_rate)
@inline parallel_state_functions(::AbstractCapPlasticity) = compute_stress, compute_pressure, compute_lambda, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate

@inline function compute_strain_rate(r::AbstractCapPlasticity; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    F = compute_F(r, τ, P)
    return ε_pl / 2 * (F > -1.0e-8)
end
@inline function compute_volumetric_strain_rate(r::AbstractCapPlasticity; τ = 0, λ = 0, P = 0, kwargs...)
    θ_pl = compute_volumetric_plastic_strain_rate(r; τ_pl = τ, λ = λ, P_pl = P, θ = 0, kwargs...)
    F = compute_F(r, τ, P)
    return θ_pl * (F > -1.0e-8)
end

# Consistency (Kuhn–Tucker) residual of the cap.
#
# Above yield this is the Duvaut–Lions equation `F = λ·η_vp`, so an
# unregularized cap (`η_vp = 0`) satisfies the rate-independent condition
# `F = 0` exactly at convergence. Below yield the `-F` term is switched off and
# the equation reduces to `λ·η_vp = 0`, which pins the multiplier at zero.
#
# The previous form added a bare `+ λ` outside the switch. That is
# dimensionally wrong — `F` and `λ·η_vp` are stresses while `λ` is a rate, so
# the equation changed its root under a change of time unit — and it also
# biased the yielding branch, whose converged state solved `F = λ·(η_vp + 1)`
# rather than `F = λ·η_vp`. Its only real job was to keep the below-yield
# branch non-degenerate when `η_vp = 0`; `_lambda_below_yield` does that
# without introducing a unit-dependent term.
@inline function compute_lambda(r::AbstractCapPlasticity; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    yielding = F > -1.0e-8
    return -F * yielding + λ * _lambda_below_yield(r, yielding)
end

# Stress-per-unit-rate factor multiplying `λ`. This is `η_vp` wherever the model
# is regularized. An unregularized cap has no viscosity at all, so below yield
# `λ·η_vp = 0` would hold for every `λ` and leave the multiplier undetermined;
# there the factor falls back to `oneunit`, which restores the old form's
# ability to pin `λ = 0` on a branch where `λ` is zero at the solution anyway,
# so no unit-dependent quantity enters the converged state. Above yield the
# factor is always `η_vp`, which is what makes `F = 0` exact.
@inline _lambda_below_yield(r::AbstractCapPlasticity, yielding) =
    (iszero(r.η_vp) && !yielding) ? oneunit(r.η_vp) : r.η_vp

@inline compute_stress(r::AbstractCapPlasticity; τ_pl = 0, kwargs...) = τ_pl
@inline compute_pressure(r::AbstractCapPlasticity; P_pl = 0, kwargs...) = P_pl

# The flow rule is the gradient of the potential Q, taken by automatic
# differentiation so a subtype needs only to define Q itself.
@inline function compute_plastic_strain_rate(r::AbstractCapPlasticity; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ * ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) - 0 * ε
end

@inline function compute_volumetric_plastic_strain_rate(r::AbstractCapPlasticity; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return -λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) - 0 * θ
end
