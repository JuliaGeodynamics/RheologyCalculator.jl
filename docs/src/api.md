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
