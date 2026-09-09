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

    J = jacobian(c, x, vars, others)
    dRdε = ForwardDiff.derivative(e -> compute_residual(c, x, merge(vars, (; ε = e)), others), εII)
    dxdε = backsolve(J, dRdε)

    return dxdε[stress_index(c)]
end

tangent(c::AbstractCompositeModel, sol::RCSolution, vars0, others) = tangent(c, sol.x, vars0, others)

"""
    tangent_block(c, x, vars, others)

Return the scalar deviatoric/volumetric consistent tangent block

```text
[ dτ/dε  dτ/dθ ]
[ dP/dε  dP/dθ ]
```

using the package's implicit local solve. Missing pressure or volumetric
unknowns produce zero rows or columns. Plastic coupling is included whenever
the model contributes the corresponding residual equations.
"""
function tangent_block(c::AbstractCompositeModel, x::SVector, vars0, others)
    vars0.ε isa Number || throw(ArgumentError("tangent_block requires scalar `vars.ε`"))
    vars = merge(vars0, (; ε = second_invariant_value(vars0.ε .+ _direct_leaf_elastic_correction(c, vars0.ε, others))))
    J = jacobian(c, x, vars, others)
    ε = vars.ε
    dRdε = ForwardDiff.derivative(e -> compute_residual(c, x, merge(vars, (; ε = e)), others), ε)
    dxdε = backsolve(J, dRdε)
    dxdθ = if hasproperty(vars, :θ)
        θ = vars.θ
        dRdθ = ForwardDiff.derivative(t -> compute_residual(c, x, merge(vars, (; θ = t)), others), θ)
        backsolve(J, dRdθ)
    else
        zero(dxdε)
    end
    τi = stress_index(c)
    Pi = findfirst(==(:P), x_keys(c))
    τ_row = (dxdε[τi], dxdθ[τi])
    P_row = Pi === nothing ? (zero(dxdε[τi]), zero(dxdε[τi])) : (dxdε[Pi], dxdθ[Pi])
    return SMatrix{2, 2}(τ_row[1], P_row[1], τ_row[2], P_row[2])
end

tangent_block(c::AbstractCompositeModel, sol::RCSolution, vars, others) =
    tangent_block(c, sol.x, vars, others)

"""
    tangent_tensor(c, x, vars, others)

Return the isotropic deviatoric Voigt tangent `∂τᵢ/∂εⱼ` at a converged local
solution. The supported tensor layouts are the package conventions `(xx, yy,
xy)` and `(xx, yy, zz, yz, xz, xy)`. The result is a fixed-size `SMatrix`.

This first form covers deviatoric invariant response only. Volumetric and
plastic coupling blocks are intentionally outside this API.
"""
function tangent_tensor(c::AbstractCompositeModel, x::SVector, vars0, others)
    ε_corr = _direct_leaf_elastic_correction(c, vars0.ε, others)
    ε = vars0.ε .+ ε_corr
    ε isa NTuple || throw(ArgumentError("tangent_tensor requires a Voigt strain-rate tuple"))
    N = length(ε)
    N ∈ (3, 6) || throw(ArgumentError("tangent_tensor supports 2D or 3D Voigt tuples"))
    εII = second_invariant_value(ε)
    iszero(εII) && throw(ArgumentError("tangent_tensor is undefined at zero strain rate"))
    τ = x[stress_index(c)]
    dτdε = tangent(c, x, vars0, others)
    εv = SVector{N}(ε)
    f = y -> begin
        e = second_invariant_value(Tuple(y))
        τy = τ + dτdε * (e - εII)
        (τy / e) .* y
    end
    return SMatrix{N, N}(ForwardDiff.jacobian(f, εv))
end

tangent_tensor(c::AbstractCompositeModel, sol::RCSolution, vars, others) =
    tangent_tensor(c, sol.x, vars, others)

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
