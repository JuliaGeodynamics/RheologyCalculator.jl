# normalisation factors for the local x vector

"""
    normalisation_x(c::AbstractCompositeModel, char_τ=1.0, char_ε=1.0)
    normalisation_x(eqs, char_τ, char_ε)

Return an `SVector` of normalization factors matching the solver vector layout
for composite model `c` or equation tuple `eqs`.

Stress-like unknowns (`τ`, `P`, `λ`) use `char_τ`; strain-rate-like unknowns
(`ε`, `θ`, plastic strain rates) use `char_ε`.

A characteristic scale of zero is replaced by one. A zero scale would make
every row it applies to unmeasurable: [`mynorm`](@ref) divides by these
factors, so a zero factor drops its row from the residual norm entirely. At a
point that is not being deformed, `char_ε = εII + |θ|` is legitimately zero,
and every strain-rate row would then be dropped — for a composite without a
plasticity element that is *all* of them, and [`solve`](@ref) would report a
residual of exactly `0.0` and stop after one iteration no matter how large the
true residual was. Since these factors only set the units in which the residual
is measured, substituting one for a vanishing scale keeps every row counted and
changes nothing for a state whose scales are nonzero.
"""
function normalisation_x(c::AbstractCompositeModel, char_τ = 1.0, char_ε = 1.0)
    eqs = generate_equations(c)
    x0 = normalisation_x(eqs, char_τ, char_ε)
    return SA[x0...]
end

@inline normalisation_x(eqs::NTuple{N, CompositeEquation}, char_τ, char_ε) where {N} =
    maptuple(eq -> _normalize_x_value(eq.fn, _nonzero_scale(char_τ), _nonzero_scale(char_ε)), eqs)

# Replace a vanishing (or non-finite) characteristic scale by one; see above.
@inline _nonzero_scale(s) = (iszero(s) || !isfinite(s)) ? oneunit(s) : abs(s)

for fn in (:compute_stress, :compute_pressure, :compute_lambda, :compute_lambda_parallel)
    @eval _normalize_x_value(::typeof($fn), char_stress, char_strainrate) = char_stress
end

for fn in (:compute_strain_rate, :compute_volumetric_strain_rate, :compute_plastic_strain_rate, :compute_volumetric_plastic_strain_rate)
    @eval _normalize_x_value(::typeof($fn), char_stress, char_strainrate) = char_strainrate
end

"""
    correct_xnorm(x, xnorm)

Return a normalization vector compatible with `x`. If `xnorm === nothing`, a
static vector of ones is used.
"""
@inline correct_xnorm(::SVector{N}, xnorm::SVector{N}) where {N} = xnorm
@inline function correct_xnorm(::SVector{N}, xnorm::SVector{M}) where {N, M}
    throw(DimensionMismatch("xnorm0 has $M entries but the solver vector has $N"))
end
@inline correct_xnorm(::SVector{N, T}, ::Nothing) where {N, T} = @SVector ones(T, N)
