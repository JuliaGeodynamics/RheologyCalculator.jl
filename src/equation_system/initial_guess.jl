# Initial guess for the local solution vector

"""
    x0 = initial_guess_x(c, vars, args, others)

Compute the initial guess `x0` for the non-linear solver given the composite rheology
model `c`. Internally generates the equation set from `c` and delegates to the
`NTuple{N, CompositeEquation}` method.

# Arguments
- `c::AbstractCompositeModel`: composite rheology model (e.g. `SeriesModel`, `ParallelModel`).
- `vars::NamedTuple`: kinematic input variables held constant during the solve
  (e.g. `(; ε = εᵢⱼ, θ = θ)` for deviatoric and volumetric strain-rate components).
- `args::NamedTuple`: current values of the differentiable unknowns that the solver
  iterates on (e.g. `(; τ = τ, P = P)`).
- `others::NamedTuple`: non-differentiable auxiliary variables required by the state
  functions, such as the time step `dt`, previous elastic stresses `τ0`, previous
  elastic pressure `P0`, or grain size `d`.

# Returns
- `x0::SVector`: static vector of initial guesses, one entry per equation in `c`.
"""
function initial_guess_x(c, vars, args, others)
    eqs = generate_equations(c)
    x0 = initial_guess_x(eqs, vars, args, others)
    return SA[promote(x0...)...]
end

"""
    x0 = initial_guess_x(eqs, vars, args, others)

Compute the initial guess for the local solution vector `x` given a tuple `eqs` of
`CompositeEquation`s. Each component of `x0` is estimated independently by calling
`estimate_initial_value` on the corresponding equation.

# Arguments
- `eqs::NTuple{N, CompositeEquation}`: equations generated from an `AbstractCompositeModel`,
  one per unknown in the solution vector.
- `vars::NamedTuple`: kinematic input variables held constant during the solve
  (e.g. `(; ε = εᵢⱼ, θ = θ)`).
- `args::NamedTuple`: current values of the differentiable unknowns (e.g. `(; τ = τ, P = P)`).
- `others::NamedTuple`: non-differentiable auxiliary variables (e.g. `dt`, `τ0`, `P0`, `d`).

# Returns
- `x0::NTuple{N}`: tuple of scalar initial guesses, one per equation.
"""
@generated function initial_guess_x(eqs::NTuple{N, CompositeEquation}, vars, args, others) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> estimate_initial_value(eqs[i], vars, args, others)
    end
end

"""
    x_keys = x_keys(c::AbstractCompositeModel)

Return the keys of the local solution vector `x` for the composite model `c`.
Keys correspond to the differentiable unknowns and may be repeated when multiple
equations share the same physical unknown (e.g. `τ` appearing in both deviatoric
and volumetric equations).

# Arguments
- `c::AbstractCompositeModel`: composite rheology model.

# Returns
- Flattened `NTuple` of `Symbol`s (e.g. `(:τ, :P, :ε, ...)`), one per equation.
"""
x_keys(c::AbstractCompositeModel) = x_keys(generate_equations(c))

"""
    x_keys = x_keys(eqs::NTuple{N, CompositeEquation})

Return a flattened tuple of `Symbol`s corresponding to the differentiable keyword argument
keys of each equation in `eqs`. Keys may be repeated when multiple equations share the
same unknown.

# Arguments
- `eqs::NTuple{N, CompositeEquation}`: equations generated from an `AbstractCompositeModel`.

# Returns
- Flattened `NTuple` of `Symbol`s (e.g. `(:τ, :P, :ε, ...)`), one per equation.
"""
@generated function x_keys(eqs::NTuple{N, CompositeEquation}) where {N}
    return quote
        @inline
        k = Base.@ntuple $N i -> begin
            keys(differentiable_kwargs(eqs[i].fn))
        end
        superflatten(k)
    end
end

"""
    estimate_initial_value(eq::CompositeEquation, vars, args, others)

Dispatch the initial-value estimation for a single equation `eq` to the appropriate
method based on the equation's kernel function (`eq.fn`):
- `compute_strain_rate` / `compute_volumetric_strain_rate` → harmonic mean over rheology components.
- `compute_stress` / `compute_pressure` → arithmetic sum over rheology components.
- any other function → returns `0`.

# Arguments
- `eq::CompositeEquation`: a single equation from the composite model, carrying the kernel
  function `eq.fn`, the tuple of rheology components `eq.rheology`, and the element
  indices `eq.el_number`.
- `vars::NamedTuple`: kinematic input variables (e.g. `(; ε = εᵢⱼ, θ = θ)`).
- `args::NamedTuple`: current differentiable unknowns (e.g. `(; τ = τ, P = P)`).
- `others::NamedTuple`: non-differentiable auxiliary variables (e.g. `dt`, `τ0`, `P0`).

# Returns
- `Float64`: scalar initial-guess value for the unknown of `eq`.
"""
estimate_initial_value(eq::CompositeEquation, vars, args, others) = _estimate_initial_value(eq.fn, eq, vars, args, others)
# Fallback: unknown equation type → use 0 as the initial guess.
@inline _estimate_initial_value(::F, eq, vars, args, others) where {F} = 0
# Strain-rate-like unknowns use a harmonic-mean estimate across the element rheologies.
@inline _estimate_initial_value(::typeof(compute_volumetric_strain_rate), eq, vars, args, others) = _estimate_initial_value_harm(eq.fn, eq.rheology, eq.el_number, vars, args, others)
@inline _estimate_initial_value(::typeof(compute_strain_rate), eq, vars, args, others) = _estimate_initial_value_harm(eq.fn, eq.rheology, eq.el_number, vars, args, others)
# Stress-like unknowns use an arithmetic-sum estimate across the element rheologies.
@inline _estimate_initial_value(::typeof(compute_pressure), eq, vars, args, others) = _estimate_initial_value_arith(eq.fn, eq.rheology, eq.el_number, vars, args, others)

# The unknown of a `compute_stress` equation is the strain rate of a parallel
# branch. The arithmetic estimate sums `compute_strain_rate` over the branch
# elements at the supplied `args.τ`, which vanishes for the usual τ = 0 seed.
# Zero is not a usable seed: it is a singular point of the parallel effective
# viscosity of a power law (η = τ/(2ε̇) ∝ τ^(1-n) → ∞ as τ → 0), so the first
# Newton step is NaN. The prescribed strain-rate invariant is finite and is the
# physical order of magnitude of a branch strain rate.
@inline function _estimate_initial_value(::typeof(compute_stress), eq, vars, args, others)
    est = _estimate_initial_value_arith(eq.fn, eq.rheology, eq.el_number, vars, args, others)
    return iszero(est) ? _prescribed_strain_rate(vars) : est
end

# `vars.ε` may be a tensor or an invariant; both reduce here. A composite with no
# prescribed deviatoric strain rate keeps the previous zero seed.
@inline _prescribed_strain_rate(vars::NamedTuple) =
    hasfield(typeof(vars), :ε) ? second_invariant_value(vars.ε) : 0

# Base cases for empty rheology tuples.
@inline _estimate_initial_value_harm(fn, rheology::Tuple{}, el_number, vars, args, others) = 1
@inline _estimate_initial_value_arith(fn, rheology::Tuple{}, el_number, vars, args, others) = 1

"""
    _estimate_initial_value_harm(fn, rheology, el_number, vars, args, others)
    _estimate_initial_value_arith(fn, rheology, el_number, vars, args, others)

Estimate an initial value for one unknown by evaluating `counterpart(fn)`
independently on each rheology component and combining the results:

| Entry point | Unknown | Combination |
|---|---|---|
| `_estimate_initial_value_harm`  | strain-rate-like | harmonic mean |
| `_estimate_initial_value_arith` | stress-like      | arithmetic sum |

Elements in series share a stress and add their strain rates, so a strain-rate
unknown combines harmonically; elements in parallel share a strain rate and add
their stresses, so a stress unknown combines arithmetically.

For each component `i`, history-dependent kwargs (e.g. `τ0`, `P0`) are extracted
from `others` using `el_number[i]` and merged with `args` and `vars`. The
combined `NamedTuple` is reduced to scalar invariants by `tensor2invariant`
before `counterpart(fn)` is called.

# Arguments
- `fn`: the kernel state function whose counterpart supplies the estimate
  (e.g. `compute_strain_rate` → `compute_stress`).
- `rheology::NTuple{N, AbstractRheology}`: the components of the model.
- `el_number`: one integer index per component, selecting its element-local
  history entries from `others`.
- `vars::NamedTuple`: kinematic inputs (e.g. `(; ε = εᵢⱼ, θ = θ)`).
- `args::NamedTuple`: the current differentiable unknowns (e.g. `(; τ, P)`).
- `others::NamedTuple`: nondifferentiable auxiliaries (e.g. `dt`, `τ0`, `P0`, `d`).

# Returns
`1` for an empty rheology tuple, and for the harmonic estimate also when every
component value is zero. Otherwise the combined estimate.
"""
@inline _estimate_initial_value_harm(fn, rheology::NTuple{N, AbstractRheology}, el_number, vars, args, others) where {N} =
    safe_inv_one(_accumulate_initial_value(safe_inv, fn, rheology, el_number, vars, args, others))

@inline _estimate_initial_value_arith(fn, rheology::NTuple{N, AbstractRheology}, el_number, vars, args, others) where {N} =
    _accumulate_initial_value(identity, fn, rheology, el_number, vars, args, others)

# Sum `weight(counterpart(fn)(component))` over the components. Unrolled so each
# component's history extraction and state-function call specialize on its own
# concrete type.
@generated function _accumulate_initial_value(weight::W, fn, rheology::NTuple{N, AbstractRheology}, el_number, vars, args, others) where {W, N}
    return quote
        @inline
        sum_vals = 0.0
        Base.@nexprs $N i -> begin
            keys_hist = history_kwargs(rheology[i])
            args_local = extract_local_kwargs(others, keys_hist, el_number[i])
            args_combined = merge(args, args_local, vars)
            fn_c = counterpart(fn)
            args_invariant = tensor2invariant(args_combined)
            sum_vals += weight(fn_c(rheology[i], args_invariant))
        end
        return sum_vals
    end
end

"""
    tensor2invariant(x)

Reduce tensor components to their scalar second invariant, dispatching on the type of `x`:

| Type of `x`                          | Behaviour                                                          |
|:-------------------------------------|:-------------------------------------------------------------------|
| `NTuple{N, Real}`                    | Returns `second_invariant(x...)` — scalar invariant of a flat component tuple (e.g. `(τxx, τyy, τxy)`). |
| `NTuple{N1, NTuple{N2, Real}}`       | Applies `tensor2invariant` element-wise; returns an `NTuple{N1}` of invariants, one per inner tuple. |
| `Tuple{}`                            | Returned unchanged.                                                 |
| `Number`                             | Returned unchanged; already a scalar.                              |
| `NamedTuple`                         | Applies `tensor2invariant` to every field value and reconstructs a `NamedTuple` with the same keys. |

# Arguments
- `x`: tensor data in one of the supported forms above.

# Returns
- A scalar, tuple of scalars, or `NamedTuple` of scalars depending on the input type.
"""
@inline tensor2invariant(A::Tuple{}) = A
@inline tensor2invariant(A::NTuple{N, Real}) where N = second_invariant(A...)
@generated function tensor2invariant(A::NTuple{N1, NTuple{N2, Real}}) where {N1, N2} 
    quote 
        @inline
        Base.@ntuple $N1 i -> tensor2invariant(A[i])
    end
end
@inline tensor2invariant(a::Number) = a

function tensor2invariant(x::NamedTuple)
    k = keys(x)
    v = values(x)
    invariants = ntuple(Val(length(k))) do i
        @inline 
        tensor2invariant(v[i])
    end
    (; zip(k, invariants)...)
end
