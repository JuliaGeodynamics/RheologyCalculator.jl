# AGENT.md — RheologyCalculator.jl

Orientation for an agent working in this repository. It describes what the
package is, how the pieces fit together, and the conventions a change must
respect.

## What the package does

RheologyCalculator.jl builds a **local** (single material point) constitutive
model out of small viscous, elastic, and plastic elements composed in series
and parallel, converts that composition into a nonlinear residual system, and
solves it with a Newton iteration. It is the rheology kernel for a host code
(typically a Stokes solver), which calls it once per integration point per
nonlinear iteration.

Everything is scalar-invariant based: the unknowns are second invariants
(`τII`, `εII`, …) plus pressure/volumetric quantities, and tensor input is
reduced to invariants before the solve. There is no field, no mesh, no time
loop here — history enters only through `τ0`, `P0`, and `dt`.

## Module layout

Two modules, one nested inside the other:

| Module | Role |
| --- | --- |
| `RheologyCalculator` | The engine: composition containers, equation generation, the Newton solver, the elastic-correction machinery, and the abstract state-function interface. Knows nothing about specific material laws. |
| `RheologyCalculator.RheologyModels` | The catalogue: `LinearViscosity`, `PowerLawViscosity`, `Elasticity`, `DruckerPrager`, creep laws, Cam-Clay variants, etc. Concrete elements only. |

Users normally load both. `RheologyModels` exports the basic elements;
advanced models (`DruckerPragerCap`, `Hyperbolic`, `ModCamClay`, `Golchin`,
`RateStateFriction`) are deliberately unexported and reached by prefix.

Directory map:

```
src/
  RheologyCalculator.jl         module definition, include order, exports
  RheologyModels.jl             submodule; includes every concrete element
  core/
    rheology_types.jl           AbstractRheology/Viscosity/Elasticity/Plasticity;
                                series/parallel_state_functions plumbing
    state_functions.jl          the 12 state-function names + fallbacks + docs
    composite.jl                SeriesModel, ParallelModel, leaf/branch sorting
    kwargs.jl                   differentiable_kwargs, history_kwargs, split_args
    others.jl                   superflatten, isvolumetric, safe_inv
  equation_system/
    equations.jl                CompositeEquation, generate_equations,
                                compute_residual  (the heart of the package)
    initial_guess.jl            initial_guess_x, x_keys, estimate_initial_value
    normalize_x.jl              normalisation_x, correct_xnorm
    solver.jl                   solve, bt_line_search, mynorm, backsolve
  post_processing/
    strain_rate_correction.jl   effective_strain_rate_correction + KV/generalized
                                Maxwell machinery, iselastic, count_elastic
    post_calculations.jl        compute_stress_elastic, compute_pressure_elastic
  rheology/                     concrete elements, grouped viscous/elastic/plastic
  display/print_rheology.jl     pretty-printing of composites
  utils/tensor_helpers.jl       2D/3D invariants and tensor <-> invariant helpers
ext/                            SparseConnectivityTracer extension
test/                           auto-discovered `test_*.jl` files
docs/                           Documenter + DocumenterVitepress site
examples/, prototypes/, 2D/     scripts, not part of the package or CI
```

## The central abstraction

### Elements

A concrete element subtypes `AbstractViscosity`, `AbstractElasticity`, or
`AbstractPlasticity`, and declares:

- `series_state_functions(r)` — which state functions it contributes when it
  sits in a `SeriesModel` (a viscous element: `(compute_strain_rate,)`).
- `parallel_state_functions(r)` — same for a `ParallelModel` (a viscous
  element: `(compute_stress,)`).
- the state functions themselves, as keyword methods
  (`compute_strain_rate(r::LinearViscosity; τ = 0, kwargs...) = τ/(2r.η)`).

`src/core/state_functions.jl` defines the twelve state-function names, a
zero-returning fallback for every `AbstractRheology`, and the
`fn(r, nt::NamedTuple)` wrapper that generated residual code calls. See
[linear_viscosity.jl](src/rheology/viscous/Newtonian/linear_viscosity.jl) for
the minimal template.

### Composition

`SeriesModel(...)` and `ParallelModel(...)` sort their arguments into two
fields at construction:

- `leafs` — direct `AbstractRheology` elements;
- `branches` — nested composites of the *opposite* kind (a `SeriesModel`
  stores `ParallelModel`s, and vice versa).

The sorting is done by `series_leafs`/`series_branches` and their parallel
twins in [composite.jl](src/core/composite.jl). Both fields are type
parameters, so the whole topology is known at compile time.

Series composition adds strain rates at common stress; parallel composition
adds stresses at common strain rate. Each parallel branch inside a series model
therefore introduces one extra unknown: that branch's strain rate.

### Equations

`generate_equations(c)` walks the composite and returns an `NTuple` of
`CompositeEquation`s, one per unknown. Each carries `parent`, `child`, `self`
(its own index), the state function `fn`, its rheology tuple, `ind_input`, and
per-element numbering used to index history tuples.

Two invariants hold throughout residual assembly and must not be broken:

- **equation position in the tuple equals `.self`** — `add_child` and
  `subtract_parent` index `x` by these numbers directly;
- **`length(generate_equations(c)) == length(x)`**, asserted in
  `compute_residual`.

`compute_residual(c, x, vars, others)` then does, in order:

1. `evaluate_state_functions` — each equation's own contribution;
2. `add_children` — add the unknowns of child equations;
3. `subtract_parent` — subtract the parent unknown, or, for a global equation,
   the prescribed input from `vars`;
4. (`SeriesModel` only) `subtract_elastic_correction` — the implicit backstress
   term for elastic elements inside parallel branches.

The unknown associated with a state function is given by
`differentiable_kwargs`: `compute_strain_rate → :τ`, `compute_stress → :ε`,
`compute_volumetric_strain_rate → :P`, `compute_pressure → :θ`,
`compute_lambda → :λ`, plastic variants → `:τ_pl`, `:P_pl`. `x_keys(c)` returns
this layout for a whole composite; `normalisation_x(c, char_τ, char_ε)` returns
matching scaling factors.

### Argument channels

Three named tuples flow through everything, and the distinction matters:

- `vars` — prescribed, **differentiable** inputs: `ε`, `θ`. Tensor-valued at
  the `solve` boundary, scalar invariants inside.
- `args` — current values of the unknowns (`τ`, `P`, …), assembled from `x`.
- `others` — **nondifferentiable** auxiliaries: `dt`, `τ0`, `P0`, grain size
  `d`. History entries are *tuples with one entry per element of that kind*
  (`τ0 = (τ0_1, τ0_2)` for two elastic elements) and are indexed by the
  equation's element numbering, not by tensor component.

### The elastic correction

Elasticity is handled as a backstress moved to the left-hand side rather than
as an extra unknown:

```
ε_eff = ε + τ0 / (2η),   η = G·dt
```

Two disjoint paths compute it, and double-counting is the standing hazard:

- **direct elastic leafs of the outer `SeriesModel`** — corrected *once,
  before* the Newton loop, in `solve` via
  `_direct_leaf_elastic_correction`, using full tensor arithmetic so
  non-coaxial `ε`/`τ0` pairs come out right;
- **elastic elements inside `ParallelModel` branches** (Kelvin-Voigt /
  generalized Maxwell) — corrected *implicitly, inside* `compute_residual` via
  `subtract_elastic_correction`, because the correction depends on the branch
  strain rates, which are unknowns.

`_kv_corrections` / `_η_KV` / `_weighted_backstress` implement the branch case.
They are derived for at most one level of Series/Parallel alternation;
`_assert_kv_nesting_supported` raises at specialization time if an elastic
element sits deeper than that, rather than silently omitting its contribution.

### The solver

[solver.jl](src/equation_system/solver.jl) is small and deliberately so:
Newton with a backtracking line search, `SMatrix` Jacobian from ForwardDiff, a
normalized L1 residual norm, and a direct static solve.

```julia
x = initial_guess_x(c, vars, args, others)
x = solve(c, x, vars, others; xnorm0 = normalisation_x(c, char_τ, char_ε))
```

`solve` corrects `vars.ε` to `εII` up front (see above), then iterates. It
currently **returns `x` unconditionally** — there is no convergence flag and no
error on non-convergence.

## Performance contract

This code runs per integration point, so it is written to be allocation-free
and type-stable, and that is enforced by tests:

- **`@generated` functions everywhere** — `generate_equations`,
  `evaluate_state_functions`, `_kv_corrections`, `flatten_repeated_functions`,
  the initial-guess estimators. Topology, element counts, and `τ0` indices are
  resolved at specialization time and emitted as literals.
- **`SVector`/`SMatrix` only** in the solver hot path; the Jacobian is a static
  matrix and `\` is a static solve.
- **`test_allocations.jl`** asserts zero allocations for `generate_equations`,
  `compute_residual`, and `solve`. It is skipped under coverage and run as a
  separate uninstrumented CI job (`--allocations-only`).
- **`test_type_stability.jl`** uses JET.

A change that introduces a runtime branch on topology, a `Vector`, a closure
capturing a non-constant, or a type-unstable `merge` will show up in these
tests. Check them before assuming a change is free.

Note the tension with generic-array practice: types here are small
`Tuple`s/`SVector`s indexed by compile-time-known positions, so `NTuple`
indexing is idiomatic. The `@inbounds` in the `@generated` bodies exist because
the indices are literals produced by `@nexprs`/`@ntuple`.

## Type-signature pitfalls

Two that have already produced bugs, both from over-tight annotations:

- **`NTuple{N}` means `Tuple{T,T,...}` — homogeneous.** A residual tuple whose
  entries are a mix of `Dual` and `Float64` (which is exactly what happens when
  ForwardDiff differentiates with respect to something that enters only one
  equation) does **not** match `NTuple{N}`. Use `Tuple{Vararg{Any,N}}` when
  only the length matters. This is the cause of issue #41's `MethodError`.
- **Sharing one `N` across two unrelated tuples** couples lengths that count
  different things. `effective_strain_rate_correction(..., ::NTuple{N,T},
  ::NTuple{N}, ...)` ties the tensor-component count of `ε` to the
  elastic-element count of `τ0`; a purely viscous model has `τ0 = ()` and
  cannot dispatch. This is issue #40.

## Testing

```julia
using Pkg; Pkg.test("RheologyCalculator")
```

`test/runtests.jl` auto-discovers every `test_*.jl` in `test/`, so a new file
needs no registration. Under coverage it drops `test_allocations.jl`; with
`--allocations-only` it runs *only* that file.

Existing suites: `test_VE`, `test_VEP`, `test_VEVP`, `test_VEPCap`,
`test_ModCamClay`, `test_Hyperbolic` (physics/regression), `test_jacobians`,
`test_equation_graphs` (equation generation), `test_solver_convergence`,
`test_strain_rate_correction`, `test_sparsity_tracer`, `test_allocations`,
`test_type_stability`, `test_Aqua`.

CI runs Julia 1.12 on Linux/macOS/Windows x64, plus allowed-to-fail `pre`
builds. Coverage and the allocation job run only on ubuntu-x64.

For interactive work use the Julia MCP session with Revise rather than fresh
`julia` processes; reserve `Pkg.test()` for a final pre-PR run.

## Conventions

- **Fail fast.** Prefer an error over a silently-propagated `NaN` or a quietly
  clamped value. A `NaN` leaving `solve` reappears in a host code's field
  hundreds of steps later, which is the expensive failure mode.
- **Annotate only as tightly as the implementation needs** (see the pitfalls
  above).
- **Comments state what is true now** — no history, no plan references, no
  "formerly". Issue numbers are acceptable as terse pointers.
- Docstrings on every exported symbol; `checkdocs = :exports` is on in
  `docs/make.jl`.
- New material laws: subtype the right abstract type, specialize
  `series_state_functions`/`parallel_state_functions`, implement the state
  functions, add a `test_*.jl` with an analytical or published reference.
- Julia compat is `1.12` (note: newer than the usual LTS floor — the package
  uses `public`/1.12-era features via its dependencies' compat).

## Public API surface

Exported from `RheologyCalculator`: `AbstractViscosity`, `AbstractPlasticity`,
`AbstractElasticity`, `CompositeModel`, `SeriesModel`, `ParallelModel`,
`generate_equations`, `initial_guess_x`, `x_keys`, `normalisation_x`, `solve`,
`effective_strain_rate_correction`.

Exported from `RheologyModels`: `LinearViscosity`, `PowerLawViscosity`,
`Elasticity`, `BulkElasticity`, `IncompressibleElasticity`, `BulkViscosity`,
`LTPViscosity`, `DruckerPrager`, `DiffusionCreep`, `DislocationCreep`,
`compute_stress_elastic`, `compute_pressure_elastic`.

Frequently needed but *not* exported, and therefore part of the friction
downstream users report: `compute_residual`, `iselastic`,
`second_invariant_2D`, `tensor_strain_rate_2D`, and the other tensor helpers in
`utils/tensor_helpers.jl`.

## Downstream usage that shapes requirements

The package is consumed by staggered-grid Stokes solvers. That context explains
the open issues: such a caller evaluates the rheology at every node including
stagnation points (very small `εII`), composes viscosity *bounds* as
series/parallel pairs (which is what introduces `ParallelModel` branch
unknowns), and needs a consistent tangent `dτII/dεII` for its own Newton
iteration.
