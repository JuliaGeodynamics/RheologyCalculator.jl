# [Composites](@id composites)

In `RheologyCalculator.jl` we can build composite rheologies with arbitrary
configurations in series, parallel, or hybrid nested networks.

All concrete elements in this page are provided by the material catalogue:

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

The two main constructors are [`SeriesModel`](@ref) and [`ParallelModel`](@ref).
Both store direct rheology elements as `leafs` and nested composites as
`branches`; this is the structure used internally by [`generate_equations`](@ref).

In a series model, compatible strain-rate contributions are summed and compared
with the prescribed strain-rate input. In a parallel model, compatible stress
contributions are summed and compared with the parent stress. Nested composites
alternate those roles recursively, which is how Maxwell, Kelvin-Voigt, Burgers,
and more general visco-elasto-plastic networks are represented.

## Material with a configuration in series

Example of a Maxwell visco-elastic model, with a viscous damper of viscosity
$\eta=10^{20}$ and an elastic spring with $G=10$ GPa and $K=46.67$ GPa. First
we define the individual components:

```jldoctest composites
julia> viscous_damper = LinearViscosity(1e20);

julia> elastic_spring = Elasticity(10e9, 46.67e9);
```

and then we stitch them together in a `SeriesModel` object:

```jldoctest composites
julia> maxwell_material = SeriesModel(viscous_damper, elastic_spring);
```

Any number of individual rheology models can be placed inside one
`SeriesModel`.

## Material with a configuration in parallel

Now let's define the Kelvin-Voigt visco-elastic model with the same individual
elements. In this case, put them together in a `ParallelModel`:

```jldoctest composites
julia> KelvinVoigt_material = ParallelModel(viscous_damper, elastic_spring);
```

As before, a `ParallelModel` can contain any number of rheology elements.

## Hybrid series/parallel material

We can also combine any number of series and parallel rheologies. In the
following example we define the Kelvin representation of the
[Burgers viscoelastic material](https://en.wikipedia.org/wiki/Burgers_material):

```jldoctest composites
julia> parallel_element = ParallelModel(viscous_damper, elastic_spring);

julia> Burgers_material = SeriesModel(viscous_damper, elastic_spring, parallel_element);
```

The resulting equation layout can be read off with [`inspect`](@ref), which
returns one entry per entry of the solver vector:

```jldoctest composites
julia> inspect(Burgers_material)
4-element ModelInspection:
  index  var  equation                        scope   elements
      1  τ    compute_strain_rate             global  LinearViscosity 1, Elasticity 1
      2  ε    compute_stress                  branch  LinearViscosity 2, Elasticity 2
      3  P    compute_volumetric_strain_rate  global  LinearViscosity 1, Elasticity 1
      4  θ    compute_pressure                branch  LinearViscosity 2, Elasticity 2
```

The `scope` column separates the equations of the outer `SeriesModel` from those
local to the parallel branch, and `elements` numbers the elements as
`display(Burgers_material)` draws them. The bare names alone are available as
[`x_keys`](@ref), which is useful when choosing an initial `args` tuple or a
normalization vector with [`normalisation_x`](@ref):

```jldoctest composites
julia> x_keys(Burgers_material)
(:τ, :ε, :P, :θ)
```

[`generate_equations`](@ref) returns the residual equations themselves, in the
same solver-vector order.

Solving a composite is covered in [Stress-time curve of a Burger's
material](@ref); [`solve`](@ref) returns an [`RCSolution`](@ref) whose entries
line up with the table above.
