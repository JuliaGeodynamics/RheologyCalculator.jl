# Repository context

RheologyCalculator.jl is a Julia package for composing viscous, elastic, and
plastic elements into series, parallel, and nested networks, then solving the
resulting local nonlinear systems with Newton iterations. Read `README.md` for
the public API and `CONTRIBUTING.md` for contribution conventions.

## Code map

- `src/RheologyCalculator.jl`: core module, include order, and public exports.
- `src/core/`: abstract element types, composite containers, state-function
  interfaces, tuple utilities, and input/history handling.
- `src/equation_system/`: equation generation, initial guesses, normalization,
  Newton solves, inspection, Jacobians, and consistent tangents.
- `src/RheologyModels.jl` and `src/rheology/`: bundled material catalogue,
  organized into viscous, elastic, and plastic laws.
- `src/post_processing/`: elastic corrections, component strain-rate partitions,
  dissipation, and shear heating.
- `src/utils/tensor_helpers.jl`: tensor conventions and helper operations.
- `src/validation.jl`: host-side model and input validation.
- `ext/`: optional SparseConnectivityTracer integration; keep this dependency
  optional through the package extension mechanism.
- `test/`: regression, analytical, convergence, allocation, inference, and API
  checks. `runtests.jl` automatically includes files named `test_*.jl`.
- `docs/src/` and `docs/make.jl`: Documenter/DocumenterVitepress documentation.
- `examples/`: runnable demonstrations with a separate Julia environment.
- `prototypes/`, `2D/`, and `rheologies/`: supplementary scripts; use the included
  files under `src/` as the canonical package implementation.

## Development commands

Run these from the repository root. `Project.toml` requires Julia 1.12; CI tests
1.12 on Linux, macOS, and Windows, plus allowed-to-fail prerelease jobs.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. -e 'using Pkg; Pkg.test(; coverage=false, test_args=["--allocations-only"])'
```

Use `Pkg.test()` to provision test-only dependencies such as Aqua and JET.
Individual test files can depend on imports and helpers in `test/runtests.jl`;
do not assume each file is standalone. Coverage runs skip allocation tests, so
run allocation checks without coverage when changing numerical hot paths.
CI also checks Julia formatting with Runic; follow the surrounding style and
avoid unrelated formatting changes.

Examples use their own environment, which points to the local package:

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/Maxwell_VE.jl
```

To prepare and build documentation:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

`docs/make.jl` also contains the CI deployment call. Documentation checks exported
API coverage with `checkdocs = :exports`; document new exports accordingly.

## Implementation conventions

- Preserve the separation between the core engine and material catalogue.
  Typical user code imports both `RheologyCalculator` and
  `RheologyCalculator.RheologyModels`. Some advanced materials intentionally
  require qualification or explicit imports.
- Add material laws by subtyping the appropriate abstract rheology type,
  specializing `series_state_functions` / `parallel_state_functions`, and
  implementing the corresponding state functions. Follow a nearby material
  implementation and register its include in `src/RheologyModels.jl`.
- Extend the core generics through explicit imports in `RheologyModels`;
  accidentally creating a same-named local function bypasses core dispatch.
- Keep solver hot paths type-stable and allocation-free. Preserve static tuples,
  StaticArrays, and `isbits` results used by GPU-compatible code. Changes to
  generated functions or dispatch need allocation and inference verification.
- Keep validation host-side. `solve` is deterministic; host-side recovery belongs
  in `solve_with_retries` rather than implicit retry logic inside ordinary solves.
- `solve` returns `RCSolution`; `sol.x` is the underlying static solution vector.
  Preserve solution indexing and reuse as input to subsequent solves.
- Preserve tensor-component Voigt conventions: `(xx, yy, xy)` in 2D and
  `(xx, yy, zz, yz, xz, xy)` in 3D. Shear entries are tensor components, not
  engineering shear. Check tangent and Voigt tests when changing tensor code.
- Preserve numerical behavior during refactors. For bug fixes, add a focused
  regression and use analytical results where available. Exercise nested and
  sibling composite branches when modifying equation ordering or indexing.
- Consult `SIMPLIFICATION_PLAN.md` when working on its refactoring phases;
  it contains phase-specific constraints and verification requirements, not an
  instruction to begin unrelated refactoring.

## Validation and change scope

For every performance experiment, use `benchmark/jacobian.jl` to capture a
baseline and candidate with the same harness and environment. Preserve raw
results in `benchmark/results/`, generate a comparison with `benchmark/report.jl`,
and update `benchmark/EXPERIMENTS.md` with improvements, regressions, validation,
and the keep/reject decision. Record unsuccessful approaches too. See
`benchmark/README.md` for commands and measurement limits.

For code changes, run the relevant regression checks and the full package suite
before handing off when feasible. Solver, equation-generation, or dispatch
changes should retain allocation, type-stability/JET, and GPU-invariant coverage.
For documentation-only changes, check referenced paths and commands; numerical
tests are unnecessary unless executable examples or behavior also change.

Keep changes focused on the requested task, update affected public documentation,
and report checks performed and any checks that could not run. Do not commit
generated manifests, coverage/allocation output, or documentation build artifacts;
these are excluded by `.gitignore`.
