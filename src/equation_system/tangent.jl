# Consistent tangent of the local solve, via the implicit function theorem.

"""
    tangent(c::AbstractCompositeModel, x::SVector, vars, others)

Return the consistent tangent `dτII/dεII` of the converged local solve of
composite model `c` at solution `x`.

`x` must be a converged solution, and `vars` and `others` the same arguments
that produced it, so that `compute_residual(c, x, vars, others) ≈ 0`. Pass
`vars` exactly as it was passed to [`solve`](@ref): the elastic correction of
the prescribed strain rate is applied here as well, so the derivative is taken
with respect to the strain-rate invariant that the Newton iteration actually
solved at.

The residual `R(x, εII) = 0` defines `x` implicitly as a function of `εII`, so

```
dx/dεII = -J⁻¹ ∂R/∂εII,    J = ∂R/∂x
```

and the stress entry of `dx/dεII` is the tangent. Both derivatives are taken by
forward-mode AD, so the result is exact rather than a finite difference of the
whole solve.

`J` is rebuilt here because it is needed at the converged `x`, which is one
Newton step beyond the last Jacobian `solve` forms.

# Example
```julia
c    = SeriesModel(LinearViscosity(1.0e3), PowerLawViscosity(4.0e-4, 3))
vars = (; ε = 1.0e-3, θ = 0.0)
x    = solve(c, initial_guess_x(c, vars, (; τ = 0.0), (;)), vars, (;))
dτdε = tangent(c, x, vars, (;))
```
"""
function tangent(c::AbstractCompositeModel, x::SVector, vars0, others)
    ε_corr = _direct_leaf_elastic_correction(c, vars0.ε, others)
    εII = second_invariant_value(vars0.ε .+ ε_corr)
    vars = merge(vars0, (; ε = εII))

    J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    dRdε = ForwardDiff.derivative(e -> compute_residual(c, x, merge(vars, (; ε = e)), others), εII)
    dxdε = backsolve(J, dRdε)

    return dxdε[stress_index(c)]
end

tangent(c::AbstractCompositeModel, sol::RCSolution, vars0, others) = tangent(c, sol.x, vars0, others)

"""
    stress_index(c::AbstractCompositeModel)

Return the position of the deviatoric stress invariant `τ` in the solver vector
of composite model `c`, as laid out by [`x_keys`](@ref).
"""
@inline stress_index(c::AbstractCompositeModel) = _stress_index(x_keys(c))

@inline function _stress_index(ks::NTuple{N, Symbol}) where {N}
    i = findfirst(==(:τ), ks)
    i === nothing && throw(
        ArgumentError(
            "composite has no deviatoric stress unknown; its solver vector holds $ks"
        )
    )
    return i
end
