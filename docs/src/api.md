# API

## Bundled Material Catalogue

Concrete material elements are defined in `RheologyCalculator.RheologyModels`.
Load them with `using RheologyCalculator.RheologyModels`.

```@docs
RheologyCalculator.RheologyModels
RheologyCalculator.RheologyModels.LinearViscosity
RheologyCalculator.RheologyModels.PowerLawViscosity
RheologyCalculator.RheologyModels.LTPViscosity
RheologyCalculator.RheologyModels.Elasticity
RheologyCalculator.RheologyModels.IncompressibleElasticity
RheologyCalculator.RheologyModels.BulkElasticity
RheologyCalculator.RheologyModels.BulkViscosity
RheologyCalculator.RheologyModels.DruckerPrager
RheologyCalculator.RheologyModels.DiffusionCreep
RheologyCalculator.RheologyModels.DislocationCreep
```

The cap and critical-state plastic models, and the rate-and-state friction law,
are not exported; reach them as
`using RheologyCalculator.RheologyModels: ModCamClay` and so on.

```@docs
RheologyCalculator.RheologyModels.DruckerPragerCap
RheologyCalculator.RheologyModels.Hyperbolic
RheologyCalculator.RheologyModels.ModCamClay
RheologyCalculator.RheologyModels.Golchin
RheologyCalculator.RheologyModels.RateStateFriction
```

## Composite Models

```@docs
SeriesModel
ParallelModel
generate_equations
x_keys
```

## Solver

```@docs
initial_guess_x
normalisation_x
solve
NonConvergenceError
tangent
stress_index
compute_residual
effective_strain_rate_correction
RheologyCalculator.second_invariant
```

## State Function Interface

Concrete rheology elements participate in generated residual equations by
extending these methods.

```@docs
RheologyCalculator.compute_strain_rate
RheologyCalculator.compute_stress
RheologyCalculator.compute_volumetric_strain_rate
RheologyCalculator.compute_pressure
RheologyCalculator.compute_lambda
RheologyCalculator.compute_lambda_parallel
RheologyCalculator.compute_plastic_strain_rate
RheologyCalculator.compute_volumetric_plastic_strain_rate
RheologyCalculator.compute_plastic_stress
RheologyCalculator.compute_viscosity
RheologyCalculator.compute_viscosity_series
RheologyCalculator.compute_viscosity_parallel
```

## Argument Helpers

```@docs
RheologyCalculator.history_kwargs
RheologyCalculator.differentiable_kwargs
RheologyCalculator.generate_args_template
RheologyCalculator.extract_local_kwargs
```

## Tensor Helpers

Deviatoric tensors are flat tuples in Voigt order: `(xx, yy, xy)` in 2D and
`(xx, yy, zz, yz, xz, xy)` in 3D. These build the `vars` and `τ0` arguments
[`solve`](@ref) expects, and advance the elastic history between time steps.

```@docs
RheologyCalculator.RheologyModels.second_invariant_2D
RheologyCalculator.RheologyModels.second_invariant_3D
RheologyCalculator.RheologyModels.tensor_strain_rate_2D
RheologyCalculator.RheologyModels.tensor_strain_rate_3D
RheologyCalculator.RheologyModels.vars_2D
RheologyCalculator.RheologyModels.vars_3D
RheologyCalculator.RheologyModels.zero_stress_tensor_2D
RheologyCalculator.RheologyModels.zero_stress_tensor_3D
RheologyCalculator.RheologyModels.stress_tensor_from_invariant_2D
RheologyCalculator.RheologyModels.stress_tensor_from_invariant_3D
RheologyCalculator.RheologyModels.elastic_stress_history_2D
RheologyCalculator.RheologyModels.elastic_stress_history_3D
```

## Post-Processing

```@docs
RheologyCalculator.compute_stress_elastic
RheologyCalculator.compute_pressure_elastic
```

## Internals

These functions are mostly useful when extending the package or debugging the
generated residual system.

```@docs
RheologyCalculator.CompositeEquation
RheologyCalculator.evaluate_state_function
RheologyCalculator.superflatten
RheologyCalculator.isvolumetric
RheologyCalculator.iselastic
```
