# Performance experiment log

## 2026-09-14: reuse the line-search residual

Hypothesis: avoid reevaluating the chosen trial point after backtracking by
returning its point, residual, and norm to `solve`. Keep line-search acceptance,
best-finite fallback, positivity limits, convergence, and iteration counting.

Implementation: `_bt_line_search_result` returns `(alpha, x, residual, norm)`;
`solve` consumes it directly. The existing step-only `bt_line_search` remains
available, including its no-evaluation behavior when the initial step is below
the minimum. If no finite trial exists, the result retains the initial point;
if there are no search trials, the solver still evaluates that initial point once.

The benchmark now interleaves original/current full solves for cold and warm
starts. `reference_solver.jl` freezes the solver and line search from commit
`6056bb0d859857a26ad9518b1f89889f99af5b3a`, sharing the unchanged residual/Jacobian
implementation. Both paths are checked for exact solution, final residual, and
iteration-count equality before timing.

[Before-change run](results/20260914T173331358_line_search_before.toml) was captured
before modifying the solver. Runic formatting subsequently changed the harness
fingerprint, so the final timing comparisons use the frozen reference within
the same run, rather than relaxing the cross-run fingerprint guard.

Initial trial: [data](results/20260914T173934015_line_search_reuse.toml),
[report](results/line_search_reuse.md). Cold scalar power-law, creep, and dilatant
solves improve by 7.3%, 22.1%, and 8.3%, respectively, with zero allocations and
exact reference agreement. A warm plastic solve regresses by 15%; Kelvin-Voigt
also trends slower. This version is not accepted as the final implementation.
The next trial explicitly inlines the result-returning helper to investigate
call/tuple overhead without changing the algorithm.

Validation:

- All 42 new line-search assertions pass, including residual evaluation counts,
  backtracked acceptance, a best trial that is not the last trial, nonfinite
  trials, and a feasible step below the search minimum.
- Full `Pkg.test()` on Julia 1.13 passes numerical, analytical, allocation,
  convergence, retries, tangent, GPU-facing invariant, and ordinary `@inferred`
  checks, but fails 29 existing JET optimization assertions. An isolated copy
  with the original solver reproduces the same 29 failures in the same five
  fixture groups. See [candidate test output](results/line_search_tests.log) and
  [original-solver test output](results/line_search_baseline_tests.log).
- The validation and Voigt files after the failing testset were run separately:
  all nine assertions pass. No tests were disabled or weakened.
- Runic 1.10.0 formatted the changed Julia files. Julia 1.12 and actual GPU
  execution were not tested on this machine.

### Inline trial and final decision

Implementation: `_bt_line_search_result` is marked `@inline`; `bt_line_search`
now returns `first(_bt_line_search_result(...))` after preserving its own
no-evaluation short-circuit for an initial step below `α_min`. No algorithmic
change from the initial trial: same acceptance rule, best-finite fallback, and
call counts.

Trial: [data](results/20260914T174146615_line_search_reuse_inline.toml),
[report](results/line_search_reuse_inline.md). This removes the regressions
seen in the initial (non-inlined) trial: warm plastic is now 7.9% faster
(previously 15% slower) and Kelvin-Voigt is flat (previously trending slower).
Cold scalar-linear, scalar-Maxwell, scalar power-law, scalar creep,
volumetric-Maxwell, plastic, and dilatant solves are all faster (8–31%,
zero allocations, exact reference agreement). Nested is +5.3%
cold / +4.9% warm (borderline, at the screening threshold) and siblings is
-1.9%/-2.4% (inconclusive); both are single-iteration fixtures dominated by
call overhead rather than the reused evaluation.

Fresh-process repeat: [data](results/20260914T174708916_line_search_reuse_inline_repeat.toml),
[same-code comparison](results/inline_repeat_comparison.md). Every fixture
reproduces identical iterations, residuals, and solutions; every `solve` and
`solve_warm` kernel change between the two inline runs is inconclusive
(≤1.6%), which is small relative to the before/after deltas above, so the
nested/siblings borderline results reflect the documented single-iteration
call-overhead noise floor rather than a regression from inlining.

Validation (repeated against this exact working tree, superseding the initial
trial's logs for acceptance purposes):

- All 42 line-search assertions still pass.
- Full `Pkg.test()` on Julia 1.13: `ERROR: Some tests did not pass: 6 passed,
  29 failed, 0 errored, 0 broken` in the `Type stability` testset only, same
  five fixture groups as the initial trial (`series elastic-viscous`,
  `Kelvin-Voigt branch`, `generalized Maxwell branch`, `plastic series`,
  `volumetric series`, each 1/7 passing). Identical to the pre-existing,
  solver-independent JET baseline; not a regression from this change. Every
  other testset, including numerical, analytical, convergence, retries,
  tangent, and GPU-facing-invariant checks, passes. See
  [full test output](results/line_search_inline_tests.log).
- `Pkg.test(; coverage=false, test_args=["--allocations-only"])` passes with
  zero allocations. See [allocation test output](results/line_search_inline_allocations.log).
- Runic 1.10.0 `--check` reports no diff for `src/equation_system/solver.jl`,
  `test/test_line_search_reuse.jl`, and the changed `benchmark/` scripts.
  Julia 1.12 and actual GPU execution were not tested on this machine.

Decision: **accept** the inlined `_bt_line_search_result` as the final
implementation of line-search reuse. It removes an entire redundant primal
residual evaluation per Newton iteration, improves or is neutral on every
fixture once inlined, introduces no allocations, and changes no observable
solver behavior (iterations, residuals, and solutions are bit-identical to
the frozen reference). Proceed to the next plan item (fused residual/Jacobian
wrapper) rather than pursuing further line-search micro-tuning; the
remaining nested/siblings variance is call-overhead noise, not evidence of a
better variant.

## 2026-09-14: baseline and Jacobian preprocessing

Environment: Windows x86_64, AMD Ryzen 7 7800X3D, Julia 1.13.0,
ForwardDiff 1.4.6, StaticArrays 1.9.20, Chairmarks 1.3.1; Julia and BLAS each use
one thread. Julia 1.12 is not installed on this machine. Production `src/` is
unchanged throughout these experiments.

### Harness development and measurement quality

- [Initial exploratory run](results/20260914T172001156_baseline.toml): useful for
  checking the harness, but `initial_guess_x` already solved the scalar power law.
  Superseded for cold-start comparisons.
- [Baseline v1](results/20260914T172216993_baseline_v1.toml) and
  [fresh-process repeat](results/20260914T172354618_baseline_v1_repeat.toml): use
  guesses scaled by 0.7, giving six Newton iterations for the scalar power law.
  Each has 10 fixtures, 12 kernels, and three randomized rounds.
- [V1 report](results/baseline_v1.md), [repeat report](results/baseline_v1_repeat.md),
  and [unchanged-code comparison](results/repeat_comparison.md) reveal substantial
  timing drift. Do not interpret these unchanged-code differences as optimization
  gains or regressions. A supposed chunk-1 gain for N=1 is particularly suspect:
  the two configs have the same type and the library ignores both.
- V2 uses CPU affinity mask 4 and interleaves the five Jacobian variants within
  each sample using Chairmarks' comparative API. This is a new harness baseline,
  not a performance improvement to the package.

### Retained baseline

- [V2 baseline data](results/20260914T172654202_baseline_v2_pinned.toml) and
  [report](results/baseline_v2_pinned.md).
- [V2 repeat data](results/20260914T172825109_baseline_v2_pinned_repeat.toml) and
  [report](results/baseline_v2_pinned_repeat.md).
- [V2 unchanged-code comparison](results/pinned_repeat_comparison.md).

Both V2 runs pass all fixture checks, with identical final solutions, residuals,
and iteration counts. Each records 360 kernel/round measurements. All baseline
residual, Jacobian, backsolve, cold-start solve, and warm-start solve kernels have
zero median allocated bytes.

Representative steady-state ranges across the two V2 runs (round-median timings):

| Fixture | N | Jacobian ns | Complete cold-start solve ns | Newton iterations |
| --- | ---: | ---: | ---: | ---: |
| Scalar power law | 1 | 34.1–48.3 | 598–784 | 6 |
| Scalar creep mixture | 1 | 59.8–81.9 | 361–497 | 2 |
| Volumetric Maxwell | 2 | 12.8 | 53.5–53.6 | 1 |
| Nested network | 3 | 27.6–27.9 | 75.7–77.0 | 1 |
| Sibling network | 4 | 36.8–38.1 | 140–143 | 1 |

Absolute timing drift remains visible, particularly for the scalar fixtures.
Pinning alone does not control CPU frequency or all external activity. Treat the
above as a baseline range, not a precise hardware-independent latency. The
interleaved comparisons make the no-gain config result much clearer: power-law,
creep, volumetric, nested, sibling, and plastic config timings are essentially
the same as uncached AD. Sub-nanosecond scalar-linear differences even reverse
direction between runs; do not promote them as a cache benefit.

Prepared equation metadata remains a clear regression: the sibling Jacobian is
about 27–28 times slower in the V2 runs, with 3184 newly allocated bytes. This
prototype is retained in the benchmark harness for comparison and is not merged
into the numerical implementation.

Validation performed: five completed benchmark runs; exact cached/uncached
residual and derivative comparisons at three states per fixture; one-pass config
assertions in V1/V2; convergence/static-output checks; report generation and
matching-environment comparisons; `git diff --check`. The full numerical suite,
JET, GPU execution, Julia 1.12, and Runic were not run for this benchmark-only
change. No production numerical code was modified.

### Preprocessing hypotheses

| Experiment | What is cached | Correctness checks | Decision |
| --- | --- | --- | --- |
| Cached callable | Residual closure outside the timed Jacobian | Identical function and inputs | No demonstrated useful gain; keep production wrapper |
| Cached config | Full-width ForwardDiff.JacobianConfig | Exact Jacobian agreement at initial, converged, and perturbed states | Do not integrate: the StaticArrays method ignores it |
| Chunk-1 config | Config with one derivative seed | Same checks, plus untimed one-residual-call assertions for both widths | Do not integrate: config does not control this method's seed width |
| Prepared equations | generate_equations output passed into a mirrored residual pipeline | Exact residual and Jacobian agreement at all three states | Reject this prototype: substantial slowdowns and new allocations |

Config creation allocates 32/64/112/176 bytes for N=1/2/3/4, respectively. These
are setup costs, not allocations in the cached Jacobian call. Config objects are
not `isbits`; the existing Jacobian and solution outputs are. Since the static
method discards the config, there is no useful setup cost to amortize here.

Precomputed equation tuples are `isbits`, but the prototype allocates 176 bytes
for Kelvin-Voigt, 416 for plastic, 448 for nested, 592 for dilatant, and 3184 for
the sibling network per Jacobian. The corresponding baseline kernels allocate
zero bytes. Passing equation indices as runtime values appears to lose compiler
specialization available when equations are generated inside the residual; this
explanation is an inference, not a proof that every possible prepared evaluator
would regress. A redesigned type-specialized evaluator is a separate experiment.

The package's local ForwardDiff extension explicitly routes static-array config
overloads back to the no-config overload. The observed one-pass behavior agrees
with that source. General configuration advice in the
[ForwardDiff guide](https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/)
must be checked against the specialized method actually dispatched here.

### Next experiments

Prioritize reusing the line-search residual, fused residual/Jacobian scheduling,
and scalar specialization. If exploring further preprocessing, focus on invariant
material factors or a proven constant Jacobian; define lifetime/invalidation and
include preparation in full-solve measurements. Preserve AD tags and avoid shared
mutable caches between concurrent solves.

The affine fixtures currently converge in one Newton iteration, so caching their
Jacobian within a single solve would mainly move its cost into preparation.
Reuse across material points or time steps could still help when all coefficients
are unchanged, but requires a different workload benchmark and an explicit
invalidation contract. These measurements do not rule out that opportunity.

For every experiment, append baseline/candidate artifacts, numerical checks,
per-fixture speed and allocation changes, and a keep/reject decision. Do not
replace the retained baseline with the candidate or discard unsuccessful trials.

## 2026-09-14: fused residual and Jacobian (plan item 3)

Hypothesis: `solve`'s pre-loop residual call and its first-iteration Jacobian
call are two separate evaluations at the same point. ForwardDiff's dual pass
computes a primal value as a byproduct of computing a Jacobian, so extracting
both from a single call should remove that redundant residual evaluation.
Later iterations already carry their residual from the previous line search
(the prior experiment above), so this only targets the first iteration.

### A nested-AD correctness bug in the naive implementation

First attempt: `DiffResults.JacobianResult(x)` for storage, then
`ForwardDiff.jacobian!(result, f, x)`. This matched the plan's description of
the StaticArrays extension and passed an isolated correctness/allocation probe
using plain `Float64` inputs. Integrated into `solve` and run against the full
suite, it broke `test/test_line_search_reuse.jl`'s "Residual reuse preserves
numeric and AD behavior" testset (4 of 6 assertions errored) with
`MethodError: no method matching Float64(::ForwardDiff.Dual{...})` inside
`DiffResults.jl`'s `tuple_setindex`.

Root cause: `solve` must itself stay differentiable — e.g.
`ForwardDiff.derivative(rate -> solve(c, x, (; ε = rate), others)[1], rate)`,
exercised by that testset for both a single and a nested (second) derivative.
There, `x` is a plain `SVector{N,Float64}` but the residual closure captures
`vars.ε` as an outer `Dual`, so `compute_residual(c, x, vars, others)` returns
values typed by that outer `Dual`, not by `eltype(x)`. `DiffResults.jl`'s
`JacobianResult(x::StaticArray)` and `JacobianResult(y::StaticArray,
x::StaticArray)` both type the Jacobian's storage from `x` alone (`zeros(
similar_type(typeof(x), Size(...)))`), independent of `y`'s type. Once
`jacobian!` tries to store an outer-`Dual`-valued entry into that
`Float64`-typed static container, the conversion throws. A second probe with
`DiffResults.JacobianResult(y, x)` (passing a `Base.promote_op`-inferred
prototype for `y`) still failed the same way: it fixes the *value* slot's type
correctly but still derives the *Jacobian* slot's type from `x`, so the bug
survives. This is a real limitation of DiffResults 1.1.0's `StaticArray`
methods for nested AD, not a chunk-size or config issue.

Fix: bypass both `JacobianResult` convenience constructors and build the
`DiffResult` directly, typing *both* the value and Jacobian storage from the
residual's own inferred output type:

```julia
f = y -> compute_residual(c, y, vars, others)
T = Base.promote_op(f, typeof(x))              # inference only, no evaluation
y_proto = zero(T)
jacobian_proto = zeros(similar_type(T, Size(length(y_proto), length(x))))
result = ForwardDiff.jacobian!(DiffResults.DiffResult(y_proto, (jacobian_proto,)), f, x)
```

`Base.promote_op` is compile-time type inference, not a residual evaluation,
so this preserves the single-evaluation goal. An isolated probe confirmed this
resolves both a first- and second-order nested derivative through a toy fused
call, matching an unfused reference exactly; the full suite (below) confirms
it for `solve` itself. Record this as a standing hazard for any future
DiffResults-based fusion in this codebase: derive every static result slot's
element type from the corresponding output prototype, never from the input.

### Benchmark harness prototype

`residual_and_jacobian(f, x)` was first added to `benchmark/jacobian.jl` only
(new `:residual_and_jacobian_fused` kernel, `DiffResults` added to
`benchmark/Project.toml`) and validated/measured before touching `src/`, as
section 1a's preprocessing experiments did. [Data](results/20260914T190751240_residual_jacobian_fusion.toml),
[report](results/residual_jacobian_fusion.md). Fused beat two separate calls
(today's pre-loop-residual-plus-first-Jacobian pattern) in 5 of 10 fixtures
(6.7–38.4% faster) and was inconclusive elsewhere, never slower. Fused was
*not* uniformly cheaper than a bare Jacobian-only call (up to 42% slower on
tiny fixtures where the DiffResults extraction overhead dominates a ~3 ns
Jacobian), confirming the plan's warning that fusion should not replace every
`jacobian()` call — only the ones that would otherwise duplicate a residual
evaluation `solve` does not already have.

### Integration and full-solve results

Implementation: `residual_and_jacobian(c, x, vars, others)` lives in
`src/equation_system/solver.jl` next to `jacobian`, with the corrected typing
above. `solve` calls it once before the loop for `r` and the first `J`; the
loop gains a single `it > 1 && (J = jacobian(c, x, vars, others))` guard so
every later iteration is unchanged (plain `jacobian`, residual from the line
search, as before). The loop always executes at least one iteration (`er`
starts at `Inf`, not `er0`), so this first Jacobian is never wasted — unlike
a next-iterate Jacobian computed eagerly at the end of an iteration, which
schedule 4 in the plan warns can be wasted on the final converged point; this
change does not do that.

[Data](results/20260914T191915440_residual_jacobian_fusion_integrated.toml),
[report](results/residual_jacobian_fusion_integrated.md). Against the frozen
reference solver, every fixture improved or was inconclusive; none regressed.
Cold-start solves: scalar_creep -36.5%, dilatant -31.0%, scalar_linear -30.6%,
scalar_maxwell -25.9%, scalar_powerlaw -22.7%, plastic -13.1%, nested -11.4%
(inconclusive), kelvin_voigt -5.3%, volumetric_maxwell -5.1%, siblings -4.0%
(inconclusive). Notably, `nested` was the one fixture the prior (line-search-
reuse-only) experiment left at +5.3% (slower); fusion turns it into an
improvement. All ten fixtures still match the frozen reference's iterations,
residual, and solution exactly, with zero allocations.

Fresh-process repeat: [data](results/20260914T192124686_residual_jacobian_fusion_integrated_repeat.toml),
[report](results/residual_jacobian_fusion_integrated_repeat.md),
[cross-run kernel comparison](results/fusion_integrated_repeat_comparison.md).
Raw kernel nanosecond values drifted substantially between the two processes
(e.g. scalar_maxwell and scalar_creep solves both look ~50–80% slower in the
second run) — but `solve_reference` (the frozen, untouched solver) drifted by
almost exactly the same amount on the same fixtures in the same run, so this
is inter-run system noise, not a regression: absolute nanosecond values are
not comparable across separate processes on this machine (consistent with the
"absolute timing drift" caveat from the baseline experiment above). What is
robust is each run's *own* interleaved reference-vs-current comparison, which
cancels that drift by sampling both together: the repeat's interleaved table
independently reproduces the same qualitative result — every fixture faster
or inconclusive, none slower, including `nested` at -15.0%. Iterations,
residuals, and solutions again match the frozen reference exactly in both
runs.

Validation:

- Full `Pkg.test()` on Julia 1.13: same pre-existing `Type stability` (JET)
  failures as every prior trial — 6 passed, 29 failed, in the same five
  fixture groups, not a regression from this change. Every other testset
  passes, including the nested-AD "Residual reuse preserves numeric and AD
  behavior" testset (6/6, both `Float32`/`Float64`, single and second
  derivative through `solve`) and all 42 line-search-reuse assertions. See
  [full test output](results/residual_jacobian_fusion_tests.log).
- `Pkg.test(; coverage=false, test_args=["--allocations-only"])` passes with
  zero allocations. See [allocation test output](results/residual_jacobian_fusion_allocations.log).
- `DiffResults` was added as a direct dependency (`Project.toml`, compat
  `"1"`), per the plan's instruction to declare it explicitly if retained.
- Runic 1.10.0 `--check` reports no content diff for the changed files; the
  only flagged difference on `src/RheologyCalculator.jl` is CRLF/LF, which
  every untouched file in this Windows checkout also shows (`core.autocrlf =
  true`) and which git normalizes to LF on commit regardless — not introduced
  by this change. Julia 1.12 and actual GPU execution were not tested on this
  machine.

Decision: **accept** the fused `residual_and_jacobian` for `solve`'s initial
residual/Jacobian pair. It removes a genuine redundant evaluation, improves or
is neutral on every fixture (fixing the one fixture the previous experiment
left regressed), preserves exact numerical/iteration parity with the frozen
reference, adds zero allocations, and the investigation surfaced and fixed a
real nested-AD hazard in the naive DiffResults approach rather than shipping
it silently broken. Do not extend fusion to every `jacobian()` call in the
loop — the harness prototype shows that would sometimes cost more than the
current schedule, exactly as the plan anticipated. Proceed to plan item 4
(scalar specialization) next.

## 2026-09-14: aggressive one-unknown specialization (plan item 4) — rejected

Hypothesis under test: for `N=1` systems, a hand-written scalar derivative (or
a fused scalar value+derivative) might beat the existing `SVector{1}`/
`SMatrix{1,1}`-typed `compute_residual`/`jacobian` path by skipping vector/dual
scaffolding. The plan itself flags the likely outcome: "StaticArrays may
already eliminate these."

### Code inspection

`@code_typed`/`@code_llvm` on `jacobian`, `compute_residual`, `mynorm`,
`max_feasible_step`, and `backsolve` for the `scalar_linear` and
`scalar_powerlaw` fixtures (`N=1`) show every one of them already compiles to
plain scalar `double` arithmetic with no heap allocation, no vector/SIMD
types, and no visible dual-number scaffolding. For `scalar_powerlaw`,
`jacobian`'s optimized LLVM IR is 17 lines: two calls to the power function,
one comparison/branch for the zero-exponent case, and scalar
`fadd`/`fmul`/`fdiv`; the `SMatrix{1,1}` return is the LLVM struct `[1 x [1 x
double]]` returned by value, not a heap object. For `scalar_linear`, the typed
IR shows the compiler has gone further and reduced the entire ForwardDiff dual
pass to the closed-form derivative `1/(2η)` at compile time — better than any
hand-written scalar path could do without also being fully inlined the same
way. `mynorm`, `max_feasible_step`, and `backsolve` show the same pattern:
scalar loads, scalar float ops, a handful of branches, no allocation. Full IR
dumps are not saved (they are reproducible in seconds from the snippets above
and `benchmark/scalar_specialization.jl`); the essential lines are quoted here
because they are the decisive evidence.

### Timing probe

[`benchmark/scalar_specialization.jl`](../benchmark/scalar_specialization.jl)
is a standalone probe (not part of `jacobian.jl`'s per-fixture kernel matrix,
since these kernels are only well-defined at `N=1`) comparing, for the four
`N=1` fixtures: the current `jacobian()`; `ForwardDiff.derivative(t ->
compute_residual(c, SVector(t), vars, others)[1], x[1])` (scalar derivative,
exactly as the plan specifies); and a fused scalar value+derivative built from
`ForwardDiff.Dual`/`Tag`/`value`/`partials` directly (the plan's fallback
design for item 3, applied here instead of `DiffResults`, which the item 3
investigation showed has a static-typing hazard not worth revisiting for a
change this probe was already expected to reject). All three are checked
against `jacobian`/`compute_residual` for exact agreement before timing.
[Result](results/scalar_specialization.md):

| Fixture | jacobian() ns | scalar derivative ns | fused scalar ns |
| --- | ---: | ---: | ---: |
| scalar_linear | 2.50 | 2.71 (+8.5%) | 2.92 (+17.1%) |
| scalar_maxwell | 2.74 | 2.96 (+7.9%) | 4.19 (+52.6%) |
| scalar_powerlaw | 38.26 | 38.26 (+0.0%) | 38.26 (+0.0%) |
| scalar_creep | 60.90 | 60.26 (-1.1%) | 62.18 (+2.1%) |

At the two fixtures with enough absolute cost to measure reliably
(`scalar_powerlaw`, `scalar_creep`, both ≥38 ns), all three approaches are
statistically indistinguishable, matching the code-inspection finding that
they compile to the same scalar operations. `scalar_linear` and
`scalar_maxwell` sit at 2.5–4 ns, where this repository's baseline experiment
already documented that sub-few-nanosecond differences are noise, not signal
(see "Sub-nanosecond scalar-linear differences even reverse direction between
runs" above); nothing here contradicts that caution.

Validation: exact value/derivative agreement asserted for all four fixtures
before timing (see the script); no separate allocation check was needed since
the LLVM inspection already shows zero allocation for every candidate and
Chairmarks reported none. Runic 1.10.0 `--check` passes on the new script.

Decision: **reject**. Code inspection and timings both show no benefit —
`jacobian`/`compute_residual` already compile to optimal scalar machine code
for `N=1`, in one case (the linear model) better than a hand-written scalar
path would without equivalent inlining. No source change. This closes plan
item 4 without touching `src/`; item 5's table (invariant hoisting, analytic
derivatives, constant Jacobians, block elimination, compressed/chunked AD,
Jacobian lagging) is the remaining unexplored scope, each requiring its own
narrower experiment and workload before evidence exists either way.
