# Component partition

A composite model solves for a handful of unknowns, but those unknowns do not say
how the deformation is partitioned between the mechanisms that produced them. A
visco-elasto-plastic point may report a single stress while the viscous, elastic,
and plastic components accommodate wildly different fractions of the imposed strain
rate, and which one dominates changes with the loading.

[`component_partition`](@ref) reports what each individual component is doing: its
deviatoric strain rate `ε` and the stress `τ` conjugate to it, grouped by
rheological category. Within a category the `ε`, `τ`, `elements`, and `indices`
tuples are index-aligned, so entry `i` of each refers to the same component, and
`indices` gives its number in the per-type numbering that `display(c)` draws.

## Series: one stress, partitioned strain rate

```jldoctest partition
julia> using RheologyCalculator, RheologyCalculator.RheologyModels

julia> c = SeriesModel(LinearViscosity(1.0e19), IncompressibleElasticity(1.0e10),
                       DruckerPrager(1.0e6, 30.0, 0.0));

julia> vars, others = (; ε = 1.0e-13), (; dt = 1.0e10, τ0 = (0.0,));

julia> sol = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), others), vars, others);

julia> component_partition(c, sol, vars, others)
ComponentPartition:
  viscous  ε = (4.33e-14,)  τ = (866000.0,)  indices = (1,)
  elastic  ε = (4.33e-15,)  τ = (866000.0,)  indices = (1,)
  plastic  ε = (5.24e-14,)  τ = (866000.0,)  indices = (1,)
```

The three components share one stress and partition the imposed strain rate between
them, which is what a series connection means. Displayed values are rounded to
three significant digits; the stored values are exact.

## Parallel: one strain rate, partitioned stress

```jldoctest partition
julia> cp = SeriesModel(ParallelModel(LinearViscosity(1.0e19), LinearViscosity(1.0e20)));

julia> vp, op = (; ε = 1.0e-13), (; dt = 1.0e10);

julia> sp = solve(cp, initial_guess_x(cp, vp, (; τ = 1.0e6), op), vp, op);

julia> component_partition(cp, sp, vp, op)
ComponentPartition:
  viscous  ε = (1.0e-13, 1.0e-13)  τ = (2.0e6, 2.0e7)  indices = (1, 2)
  elastic  ε = ()  τ = ()  indices = ()
  plastic  ε = ()  τ = ()  indices = ()
```

Both components deform at exactly the same rate, yet the stiffer one carries ten
times the stress and therefore dissipates ten times the power. This is why
dissipation cannot be read off a strain-rate field: each component must be paired
with its own conjugate stress, never a global stress with a summed local rate.

## Why is partitioning worth looking at

Sweeping the imposed strain rate over six orders of magnitude for the
visco-elasto-plastic point above shows what the stresses alone cannot.

![](https://raw.githubusercontent.com/juliageodynamics/RheologyCalculator.jl/main/docs/assets/component_partition.png)

Below yield, stresses rise linearly with the imposed rate and the viscous/elastic
split is fixed by the ratio ``\eta / (G \Delta t)``. Once the Drucker--Prager yield
stress ``C\cos\phi`` is reached, the stresses saturate and stop carrying
information: every further increment of imposed strain rate is taken up by plastic
flow, and only the partitioning shows that handover.

The figure is produced by `examples/component_partition.jl`.