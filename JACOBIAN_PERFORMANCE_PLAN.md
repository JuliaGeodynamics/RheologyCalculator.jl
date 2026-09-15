# Jacobian and Newton-solve performance exploration

## Objective

Reduce steady-state `solve` cost, starting with `SVector{1}` systems and redundant
residual evaluations. Measure each option independently before combining winners.
This is an experimental plan; no speedups have been measured and no solver changes
are part of this document.

## Current implementation and implications

- `src/equation_system/solver.jl` wraps `ForwardDiff.jacobian` in `jacobian`.
  `solve` retains a primal residual between iterations and computes a fresh
  Jacobian at the current iterate.
- Line search evaluates trial residuals but returns only the step length. The
  caller evaluates the chosen point again to check convergence.
- `backsolve` already specializes both `SVector{1}` and `SMatrix{1,1}` Jacobians.
  Scalarization must improve differentiation or surrounding work to add value.
- Residual construction already specializes on static sizes and composite types.
  Source-level tuple/equation construction does not prove runtime overhead;
  inspect optimized code before adding caches or more generated functions.
- The locally available ForwardDiff 1.4.6 StaticArrays extension seeds all N
  directions and evaluates the residual once. Its immutable DiffResults path
  extracts both value and Jacobian from that evaluation. StaticArray methods
  ignore supplied JacobianConfig objects, so changing their chunk size alone
  will not tune this path. Verify the actual resolved version in experiments;
  the repository allows ForwardDiff versions compatible with 1.4.1.
- Tangents also call `jacobian`; retain its static matrix return contract and
  accurate derivative semantics when experimenting with solver-only policies.

## 1. Establish a reproducible baseline

Create a benchmark script in the examples environment, which already includes
Chairmarks. Record Julia/dependency versions, CPU, thread count, and settings.
Warm each model/type combination; keep construction, setup, and compilation out
of steady-state measurements. Report first-call latency separately.

Fixtures should cover:

| Family | Purpose |
| --- | --- |
| One-unknown linear viscous / incompressible Maxwell | Scalar overhead floor |
| One-unknown power-law and diffusion/dislocation mixtures | Cheap versus expensive scalar derivatives |
| Two-unknown volumetric Maxwell | Small dense baseline |
| VEP/VEVP, dilatant and cap models | Coupling and branch-dependent laws |
| Nested and sibling series/parallel networks | Increasing N and structural sparsity |
| Difficult low-strain-rate cases and warm starts | Backtracking, positivity limits, and near-convergence work |

Assert actual unknown counts rather than inferring them from model names. Include
non-coaxial tensor/history fixtures to exercise elastic correction semantics.

Measure primal residual, Jacobian, residual-plus-Jacobian, backsolve, and full
solve. Report time distributions, allocations, iterations, final residual, and
failure rates. Count primal evaluations, dual evaluations, and line-search trials
in separate diagnostic runs so counters do not contaminate timings. Compare both
fixed-state kernels and complete solves from identical initial guesses/tolerances.

## 1a. Measure Jacobian preprocessing and cache lifetime

Before integration experiments, assess preprocessing with the harness in
`benchmark/`: cached residual callable, cached full-width and chunk-1
`JacobianConfig`, and precomputed equation metadata. Record setup time and
allocations separately from reuse. The experiment log in
`benchmark/EXPERIMENTS.md` is the running record of accepted and rejected ideas;
retain baseline and candidate TOML files and generated comparisons.

For any future useful cache, define its lifetime and invalidation explicitly:
topology can follow the composite type, material coefficients must follow model
values and environmental inputs, and history/time-step terms must refresh per
solve. Seeds can encode directions but must never retain old primal iterates.
Preserve AD tags and numeric types, avoid shared mutable workspaces between
threads/device lanes, and include preparation cost when measuring a full solve.
Report the number of reuses required to amortize setup when reuse is beneficial.

## 2. Reuse line-search results

Prototype an internal line-search result containing the chosen step, point,
residual, and normalized norm. Preserve the existing step-only helper if callers
depend on it. Reuse the result for convergence and the next Newton right-hand side.

Retain the residual of the best finite trial, not merely the last trial. Cover
immediate acceptance, backtracked acceptance, best-trial fallback, no finite
trial, and initial feasible steps below the minimum search step. Where no trial
was evaluated, preserve the caller's existing evaluation and failure behavior.
Keep acceptance thresholds, stagnation logic, and iteration counting unchanged.

This is a distinct experiment: it can remove an entire primal evaluation without
changing the derivative backend.

## 3. Combined residual and Jacobian wrapper

Prototype an internal interface:

```julia
residual_and_jacobian(c, x::SVector, vars, others) -> (r, J)
```

First use ForwardDiff's supported `jacobian!` API with an immutable DiffResults
container holding static value/Jacobian storage. Capture the returned container:
immutable result updates return a new value. Initialize storage without evaluating
the residual solely to discover its shape. Check mixed numeric/output types.
If retained, declare DiffResults as a direct dependency with compatibility bounds.

Compare with a small wrapper that seeds tagged ForwardDiff duals and extracts
static values/partials directly only if the supported API shows measurable
overhead. Avoid depending on unexported StaticArrays-extension internals. Preserve
nested differentiation and tag isolation; never strip outer dual layers.

Benchmark integration schedules explicitly:

1. Existing cached primal residual plus Jacobian-only evaluation.
2. Fused evaluation at the current iterate, reusing its primal result wherever
   this actually eliminates a separate call.
3. Reused line-search residual plus Jacobian-only evaluation.
4. Fused evaluation at the accepted next point, retaining the next Jacobian.

Fusion does not automatically save a residual call: the previous iteration or
line search may already supply it. Schedule 4 computes an unnecessary Jacobian
on the final converged point. Do not differentiate every rejected line-search
trial. Measure saved primal work against extra derivative work, particularly for
one-step and warm-start solves. Preserve the current initial-iteration semantics.

## 4. Aggressive one-unknown specialization

Compare the existing static Jacobian against a scalar derivative of
`t -> compute_residual(c, SVector(t), vars, others)[1]`, then against fused scalar
value/derivative evaluation. Keep the external Jacobian shape as `SMatrix{1,1}`;
an internal Newton-step helper can consume a scalar derivative directly.

Inspect optimized IR/native code for residual scaffolding, dual construction,
matrix construction, norm, and feasibility-mask overhead. StaticArrays may already
eliminate these. Add scalar residual or Newton helpers only where code inspection
and timings demonstrate a benefit. Prefer shared convergence/line-search logic
over maintaining a second complete solver.

Preserve zero/nonfinite derivative behavior, normalization semantics, elastic
correction, and feasible-step constraints. Do not assume every one-unknown model
is a linear stress equation. Avoid changing reciprocal/multiply arithmetic to
division silently: rounding differences can alter convergence near tolerance.

## 5. Additional options, ordered by evidence and scope

| Option | Experiment and decision criterion |
| --- | --- |
| Hoist invariant material work | Identify temperature/pressure factors, history terms, or equation metadata recomputed inside residuals. Cache only quantities independent of x for the duration of a solve; preserve outer AD dependence. Keep only improvements surviving compiler inspection and full-solve timing. |
| Analytic or hybrid derivatives | Start with linear viscosity/elasticity and simple scalar power laws. Assemble exact topology coefficients and differentiate unsupported material terms with AD. Include implicit elastic corrections and plastic coupling; avoid a duplicate constitutive implementation without shared value/derivative expressions. |
| Constant Jacobians for proven affine models | Use an explicit capability/trait with a conservative fallback; reuse J and possibly its factorization only within a solve where coefficients are fixed. Do not infer constancy from size or an inactive plastic branch. |
| Block elimination | Explore eliminating linear branch unknowns or separating uncoupled deviatoric/volumetric blocks. Benchmark total solve cost and conditioning; retain full-system validation for plastic coupling. |
| Compressed forward AD | For larger sparse nested networks, derive conservative residual dependencies and color columns before the numerical loop. Keep static seeds and output for small systems; establish an N/density crossover. The current tracer extension for solve deliberately returns conservative dense dependencies, so it is not already a useful residual coloring. |
| Chunked forward AD | For larger N, explicitly prototype a path supporting narrower seed widths; the existing static method ignores config chunk sizes. Balance repeated residual evaluations against dual width, compilation, and register pressure. |
| Jacobian lagging / quasi-Newton | Separate opt-in algorithm experiment after exact optimizations. Reuse derivatives until progress degrades, then refresh. Compare total time, iterations, line searches, and failures around yielding; do not use approximate Jacobians for reported consistent tangents. |

Alternative AD backends and matrix-free Newton are deferred until benchmarks show
a system size or workload where these additional dependencies and algorithms are
justified. The initial target is small, statically sized local systems.

## Validation and adoption

- Compare residuals and Jacobians to the current backend at identical states;
  use analytical fixtures and independent finite differences at smooth points.
  At yield switches, test each side and the established branch convention rather
  than treating a central difference across the switch as an exact reference.
- Check nested AD, Float32/Float64 and supported mixed types, static shapes,
  `isbits`, zero hot-path allocations, inference/JET, and GPU invariants. Device
  execution is an additional check when hardware/tooling is available; CPU
  `isbits` checks alone do not establish device compatibility.
- Run the full package suite and allocation tests without coverage. Prioritize
  Jacobian, tangent, tensor tangent, block tangent, convergence, nonnegative
  iterate, retries, analytical regression, and elastic-correction tests.
- For scheduling-only changes, require matching residuals, iterates, line-search
  decisions, and failure diagnostics. Investigate any numerical differences from
  fused primal extraction. For analytic arithmetic changes, document rounding
  differences and enforce justified numerical tolerances plus convergence checks.
- Record per-fixture before/after results and uncertainty. Retain a change only
  when its full-solve benefit exceeds timing noise without material regressions
  in other supported cases; use narrow dispatch if the benefit is localized.
  Do not trade compilation, allocations, or robustness for a microbenchmark win
  without reporting that tradeoff.

Suggested delivery order: baseline harness; line-search reuse; fused wrapper and
integration comparisons; scalar specialization; then whichever additional option
the measurements support. Keep experiments independently reviewable and combine
only demonstrated winners. Leave the public solver API unchanged initially.

## References

- [ForwardDiff advanced guide](https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/): retrieving primal results and differentiation configuration.
- [DiffResults API](https://juliadiff.org/DiffResults.jl/stable/): static result containers and immutable update semantics.
- [ForwardDiff StaticArrays extension](https://github.com/JuliaDiff/ForwardDiff.jl/blob/v1.4.6/ext/ForwardDiffStaticArraysExt.jl): the locally inspected implementation; verify against the benchmark environment's resolved version.
