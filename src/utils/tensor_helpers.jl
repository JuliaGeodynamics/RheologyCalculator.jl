# Deviatoric tensors are stored as flat tuples in Voigt order: `(xx, yy, xy)` in
# 2D and `(xx, yy, zz, yz, xz, xy)` in 3D. The `_2D` / `_3D` suffixes name the
# expected tuple length; each pair shares one dimension-agnostic body, which
# reads the length from the tuple it is given.

const εxx_pure_shear = (1.0, -1.0, 0.0)
const εxx_pure_shear_3D = (1.0, -1.0, 0.0, 0.0, 0.0, 0.0)

"""
    second_invariant_2D(ε)
    second_invariant_3D(ε)

Second invariant of a deviatoric tensor held as a flat tuple in Voigt order,
`(xx, yy, xy)` in 2D and `(xx, yy, zz, yz, xz, xy)` in 3D.

The 2D form reconstructs `zz = -xx - yy` from incompressibility, so it is the
invariant of the full three-dimensional tensor rather than of a plane-strain
slice.

```jldoctest
julia> using RheologyCalculator.RheologyModels

julia> second_invariant_2D((1.0, -1.0, 0.0))
1.0
```
"""
second_invariant_2D(ε) = second_invariant(ε...)

@doc (@doc second_invariant_2D)
second_invariant_3D(ε) = second_invariant(ε...)

"""
    tensor_strain_rate_2D(εII; direction = εxx_pure_shear)
    tensor_strain_rate_3D(εII; direction = εxx_pure_shear_3D)

Build a strain-rate tensor of second invariant `εII`, oriented along
`direction`. `direction` is a tuple in the same Voigt order as the result and
is normalized by its own second invariant, so only its orientation matters.

The default is pure shear in the `xx`-`yy` plane.
"""
tensor_strain_rate_2D(εII; direction = εxx_pure_shear) = _tensor_strain_rate(εII, direction)

@doc (@doc tensor_strain_rate_2D)
tensor_strain_rate_3D(εII; direction = εxx_pure_shear_3D) = _tensor_strain_rate(εII, direction)

@inline function _tensor_strain_rate(εII, direction)
    directionII = second_invariant(direction...)
    return @. εII * direction / directionII
end

"""
    vars_2D(εII, θ = 0.0; direction = εxx_pure_shear)
    vars_3D(εII, θ = 0.0; direction = εxx_pure_shear_3D)

Assemble the `vars` NamedTuple `solve` expects: a strain-rate tensor of second
invariant `εII` built by [`tensor_strain_rate_2D`](@ref), plus the volumetric
strain rate `θ`.
"""
vars_2D(εII, θ = 0.0; direction = εxx_pure_shear) = (; ε = tensor_strain_rate_2D(εII; direction), θ)

@doc (@doc vars_2D)
vars_3D(εII, θ = 0.0; direction = εxx_pure_shear_3D) = (; ε = tensor_strain_rate_3D(εII; direction), θ)

"""
    zero_stress_tensor_2D()
    zero_stress_tensor_3D()

A zero stress tensor of the right length for the dimension, for seeding the
`τ0` history of a time loop.
"""
zero_stress_tensor_2D() = (0.0, 0.0, 0.0)

@doc (@doc zero_stress_tensor_2D)
zero_stress_tensor_3D() = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

"""
    stress_tensor_from_invariant_2D(τII, ε_eff)
    stress_tensor_from_invariant_3D(τII, ε_eff)

Reconstruct a stress tensor of second invariant `τII` that is coaxial with the
effective strain-rate tensor `ε_eff`, via the effective viscosity
`τII / (2 ε_effII)`. A zero `ε_eff` gives a zero stress tensor rather than a
division by zero.
"""
stress_tensor_from_invariant_2D(τII, ε_eff) = _stress_tensor_from_invariant(τII, ε_eff)

@doc (@doc stress_tensor_from_invariant_2D)
stress_tensor_from_invariant_3D(τII, ε_eff) = _stress_tensor_from_invariant(τII, ε_eff)

function _stress_tensor_from_invariant(τII, ε_eff)
    εII = second_invariant(ε_eff...)
    η_eff = iszero(εII) ? zero(τII) : τII / (2 * εII)
    return Tuple(@. 2 * η_eff * ε_eff)
end

"""
    elastic_stress_history_2D(c, τII, ε, τ0, others)
    elastic_stress_history_2D(c, x::SVector, ε, τ0, others)
    elastic_stress_history_3D(c, τII, ε, τ0, others)
    elastic_stress_history_3D(c, x::SVector, ε, τ0, others)

Advance the elastic backstress history of composite `c` by one step, returning
one stress tensor per elastic element for use as the next step's `τ0`.

Given a scalar `τII`, the result is a one-element tuple. Given a solver vector
`x`, the elastic stresses are recovered with `compute_stress_elastic`, so the
result has one entry per elastic element in `c`.
"""
elastic_stress_history_2D(c, τII, ε, τ0, others) =
    _elastic_stress_history(stress_tensor_from_invariant_2D, c, τII, ε, τ0, others)

@doc (@doc elastic_stress_history_2D)
elastic_stress_history_3D(c, τII, ε, τ0, others) =
    _elastic_stress_history(stress_tensor_from_invariant_3D, c, τII, ε, τ0, others)

function _elastic_stress_history(to_tensor::F, c, τII, ε, τ0, others) where {F}
    ε_eff = ε .+ effective_strain_rate_correction(c, ε, τ0, others)
    return (to_tensor(τII, ε_eff),)
end

function _elastic_stress_history(to_tensor::F, c, x::SVector, ε, τ0, others) where {F}
    ε_eff = ε .+ effective_strain_rate_correction(c, ε, τ0, others)
    τII_elastic = compute_stress_elastic(c, x, others)
    return ntuple(i -> to_tensor(τII_elastic[i], ε_eff), length(τII_elastic))
end
