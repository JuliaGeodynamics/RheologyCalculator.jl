# Consistent tangent of the local solve, via the implicit function theorem.

"""
    tangent(c::AbstractCompositeModel, x::SVector, vars, others)

Return the consistent tangent `dτII/dεII` of the converged local solve of
composite model `c` at solution `x`, taken along the applied strain rate.

`x` must be a converged solution, and `vars` and `others` the same arguments
that produced it, so that `compute_residual(c, x, vars, others) ≈ 0`.

For a scalar `vars.ε` this is `dτII/dε`. For a Voigt tuple it is the derivative
of `τII` with respect to `s` along `ε + s ε/εII`, i.e. with the strain rate
scaled in its own direction; the backstress history is held fixed. The full
tensor response is [`tangent_tensor`](@ref).

The residual `R(x, ε) = 0` defines `x` implicitly as a function of `ε`, so

```
dx/ds = -J⁻¹ ∂R/∂s,    J = ∂R/∂x
```

and the stress entry of `dx/ds` is the tangent. Both derivatives are taken by
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
function tangent(c::AbstractCompositeModel, x::SVector, vars, others)
    J = jacobian(c, x, vars, others)
    dRdε = _residual_derivative_along_strain_rate(c, x, vars, others, vars.ε)
    dxdε = backsolve(J, dRdε)

    return dxdε[stress_index(c)]
end

tangent(c::AbstractCompositeModel, sol::RCSolution, vars, others) = tangent(c, sol.x, vars, others)

@inline _residual_derivative_along_strain_rate(c, x, vars, others, ε::Number) =
    ForwardDiff.derivative(e -> compute_residual(c, x, merge(vars, (; ε = e)), others), ε)

function _residual_derivative_along_strain_rate(c, x, vars, others, ε::NTuple)
    εII = second_invariant_value(ε)
    iszero(εII) && throw(ArgumentError("the tangent along the applied strain rate is undefined at zero strain rate"))
    direction = ε ./ εII
    return ForwardDiff.derivative(s -> compute_residual(c, x, merge(vars, (; ε = ε .+ s .* direction)), others), zero(εII))
end

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
function tangent_block(c::AbstractCompositeModel, x::SVector, vars, others)
    vars.ε isa Number || throw(ArgumentError("tangent_block requires scalar `vars.ε`"))
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

Return the deviatoric Voigt tangent `∂τᵢ/∂εⱼ` at a converged local solution.
The supported tensor layouts are the package conventions `(xx, yy, xy)` and
`(xx, yy, zz, yz, xz, xy)`. The result is a fixed-size `SMatrix`.

The stress tensor is `τ = τII n` with `n = E/EII` and `E = ε + H`, where `H` is
the history tensor of [`effective_strain_rate_correction`](@ref). Both `τII`
and `H` depend on `ε` through the solution `x`, whose derivative
`dx/dε = -J⁻¹ ∂R/∂ε` is taken with respect to every component. Without history
the result is symmetric; with a history not coaxial with `ε` it is not in
general.

This form covers the deviatoric response only. Volumetric and plastic coupling
blocks are outside this API.
"""
function tangent_tensor(c::AbstractCompositeModel, x::SVector, vars, others)
    ε = vars.ε
    ε isa NTuple || throw(ArgumentError("tangent_tensor requires a Voigt strain-rate tuple"))
    N = length(ε)
    N ∈ (3, 6) || throw(ArgumentError("tangent_tensor supports 2D or 3D Voigt tuples"))
    εv = SVector{N}(ε)
    J = jacobian(c, x, vars, others)
    dRdε = ForwardDiff.jacobian(y -> compute_residual(c, x, merge(vars, (; ε = Tuple(y))), others), εv)
    dxdε = J \ -dRdε
    eqs = generate_equations(c)
    τi = stress_index(c)
    f = y -> begin
        xl = x + dxdε * (y - εv)
        E = assembled_strain_rate(c, eqs, xl, Tuple(y), others)
        EII = second_invariant_value(E)
        iszero(EII) && throw(ArgumentError("tangent_tensor is undefined at a zero assembled strain rate ε + H"))
        (xl[τi] / EII) .* SVector{N}(E)
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
