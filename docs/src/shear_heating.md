# Shear heating

Shear heating is the conversion of irreversible mechanical work into heat during
viscous creep and plastic deformation. It is a
volumetric source term that can locally raise temperature and thereby affect the
material strength. Reversible elastic storage is not shear heating.

## Problem definition

Consider one material point in a visco-elasto-plastic rock undergoing a prescribed
deviatoric strain rate. The point contains linear viscous creep, an elastic spring,
and Drucker--Prager plasticity connected in series. At each time step we solve for
the stress and the strain-rate contribution of each mechanism. The dissipative
viscous and plastic work can then be supplied to an energy equation as a volumetric
heat source.

The following setup solves one such increment. `η` is viscosity in Pa s, `G` is
shear modulus in Pa, `ε` is the imposed second invariant of strain rate in s⁻¹, and
`dt` is the time step in s. `τ0` is the deviatoric stress from the preceding step.

```julia
using RheologyCalculator
using RheologyCalculator.RheologyModels

η = 1.0e19                      # Pa s
G = 1.0e10                      # Pa
ε = 1.0e-13                     # s⁻¹
dt = 1.0e10                     # s

c = SeriesModel(
    LinearViscosity(η),
    IncompressibleElasticity(G),
    DruckerPrager(1.0e6, 30.0, 0.0),
)
vars = (; ε)
others = (; dt, τ0 = (0.0,))
x = initial_guess_x(c, vars, (; τ = 1.0e6), others)
sol = solve(c, x, vars, others)
```

## Dissipation partition and shear heating

[`dissipation_partition`](@ref) separates irreversible **deviatoric**
dissipation into viscous and plastic contributions. Elastic contributions are
deliberately excluded because elastic work is reversible storage. Both fields have
units of W/m³.

```julia
partition = dissipation_partition(c, sol, vars, others)

partition.viscous_Φ       # deviatoric viscous dissipation
partition.plastic_Φ       # deviatoric plastic dissipation
```

[`shear_heating`](@ref) returns only the deviatoric viscous and plastic
dissipation. Its optional Taylor--Quinney factor `β` applies to the deviatoric
plastic term; for example, `shear_heating(partition; β = 0.9)` converts 90% of
that plastic dissipation and all deviatoric viscous dissipation into heat.
