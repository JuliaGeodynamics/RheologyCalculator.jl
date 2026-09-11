# Rheologies

Concrete rheology elements are regular Julia types that subtype one of the
abstract rheology supertypes. The package's bundled implementations live in
`RheologyCalculator.RheologyModels`; import that module before constructing
one of its exported elements:

```@meta
DocTestSetup = quote
    using RheologyCalculator
    using RheologyCalculator.RheologyModels
end
```

```julia
using RheologyCalculator
using RheologyCalculator.RheologyModels
```

The abstract types themselves are defined by the core module:

```@docs
RheologyCalculator.AbstractRheology
RheologyCalculator.AbstractViscosity
RheologyCalculator.AbstractElasticity
RheologyCalculator.AbstractPlasticity
RheologyCalculator.AbstractCapPlasticity
```

Each concrete element declares the state functions it contributes in series and
parallel composition. The state-function interface is documented below and in
[API](@ref).

## Element Interface

A concrete rheology is a small immutable Julia type plus methods for:

- `series_state_functions(r)`: functions used when `r` sits in a
  [`SeriesModel`](@ref).
- `parallel_state_functions(r)`: functions used when `r` sits in a
  [`ParallelModel`](@ref).
- State functions such as `compute_strain_rate`, `compute_stress`,
  `compute_pressure`, or plastic consistency functions.
- `history_kwargs(r)`, when values from `others` should be indexed per element.

For example, a deviatoric Newtonian viscosity contributes strain rate in series
and stress in parallel. Two such dampers in series carry the same stress and
their strain rates add, so the composite has the harmonic effective viscosity
``\eta_1\eta_2/(\eta_1+\eta_2)``, which is what the last line below recovers:

```jldoctest
struct MyLinearViscosity{T} <: RheologyCalculator.AbstractViscosity
    η::T
end

RheologyCalculator.series_state_functions(::MyLinearViscosity) =
    (RheologyCalculator.compute_strain_rate,)
RheologyCalculator.parallel_state_functions(::MyLinearViscosity) =
    (RheologyCalculator.compute_stress,)

RheologyCalculator.compute_strain_rate(r::MyLinearViscosity; τ = 0, kwargs...) = τ / (2 * r.η)
RheologyCalculator.compute_stress(r::MyLinearViscosity; ε = 0, kwargs...) = 2 * r.η * ε

c    = SeriesModel(MyLinearViscosity(1.0e20), MyLinearViscosity(2.0e20))
vars = (; ε = 1.0e-14)
x0   = initial_guess_x(c, vars, (; τ = 1.0), NamedTuple())
sol  = solve(c, x0, vars, NamedTuple())

sol[1] / (2 * vars.ε)

# output

6.666666666666666e19
```

The solver passes local arguments as keywords. Unknowns come from the solver
vector `x`, prescribed inputs come from `vars`, and auxiliary values come from
`others`. Element-local history fields are selected with
[`history_kwargs`](@ref RheologyCalculator.history_kwargs); for elastic elements
the default history fields are `τ0` and `P0`, and for viscous elements the
default is `d`.

## Creep Laws

### Linear Viscosity

```julia
LinearViscosity(η)
```

For scalar deviatoric quantities:

```math
\tau = 2\eta\dot\varepsilon
```

### Power-Law Viscosity

```julia
PowerLawViscosity(η, n)
```

The current implementation uses:

```math
\dot\varepsilon = \frac{\tau^n}{2\eta}
```

### Diffusion Creep

```julia
DiffusionCreep(n, r, p, A, E, V, R)
```

where `n` is the stress exponent, `r` is the water-fugacity exponent, `p` is the
grain-size exponent, `A` is the prefactor, `E` is the activation energy, `V` is
the activation volume, and `R` is the gas constant.

### Dislocation Creep

```julia
DislocationCreep(n, r, A, E, V, R)
```

where `n` is the stress exponent, `r` is the water-fugacity exponent, `A` is the
prefactor, `E` is the activation energy, `V` is the activation volume, and `R`
is the gas constant.

## Elasticity

### Compressible Elasticity

```julia
Elasticity(G, K)
```

where `G` and `K` are the shear and bulk moduli.

### Incompressible Elasticity

```julia
IncompressibleElasticity(G)
```

where `G` is the shear modulus.

### Bulk Elements

```julia
BulkElasticity(K)
BulkViscosity(χ)
```

## Plastic Failure

### Drucker-Prager

```julia
DruckerPrager(C, ϕ, ψ[, η_vp])
```

where `C` is the cohesion, `ϕ` and `ψ` are the friction and dilation angles,
and the optional `η_vp` is the Duvaut-Lions viscoplastic regularisation
viscosity, expressed in the same viscosity units as the rest of the model.
It defaults to `1.0` for backwards compatibility. Set `η_vp = 0` to recover
the unregularised yield condition.

## State Functions

Rheologies extend these methods to participate in equation generation:

- [`compute_strain_rate`](@ref RheologyCalculator.compute_strain_rate)
- [`compute_stress`](@ref RheologyCalculator.compute_stress)
- [`compute_volumetric_strain_rate`](@ref RheologyCalculator.compute_volumetric_strain_rate)
- [`compute_pressure`](@ref RheologyCalculator.compute_pressure)
- [`compute_lambda`](@ref RheologyCalculator.compute_lambda)
- [`compute_lambda_parallel`](@ref RheologyCalculator.compute_lambda_parallel)
- [`compute_plastic_strain_rate`](@ref RheologyCalculator.compute_plastic_strain_rate)
- [`compute_volumetric_plastic_strain_rate`](@ref RheologyCalculator.compute_volumetric_plastic_strain_rate)
- [`compute_plastic_stress`](@ref RheologyCalculator.compute_plastic_stress)
- [`compute_viscosity`](@ref RheologyCalculator.compute_viscosity)
- [`compute_viscosity_series`](@ref RheologyCalculator.compute_viscosity_series)
- [`compute_viscosity_parallel`](@ref RheologyCalculator.compute_viscosity_parallel)
