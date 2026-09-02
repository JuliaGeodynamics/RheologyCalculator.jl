# Simplification plan for `src/`

Target: reduce boilerplate, duplication, and stale documentation in
RheologyCalculator's 4265 lines of `src/`, without changing a single numerical
result or removing a working feature.

Baseline verified on this machine: Julia 1.12.7, package loads, `solve` on
`SeriesModel(Elasticity, LinearViscosity)` returns `[66666.66666666667, -100.0]`.

## Ground rules

1. **Behavior-preserving, except where explicitly fixing a bug.** Phase 1 changes
   results — that is its purpose, and each change there is a defect fix with a
   reproduction. Every other phase is a refactor, a deletion of unreachable code,
   or a documentation fix, and must leave numerical output bit-identical.
2. **One phase per PR.** Phases are ordered so each rests on the verification
   net established by the one before it. Do not interleave.
3. **Verification is not optional.** The package already has the two guardrails
   that make this kind of refactor safe — `test/test_allocations.jl`
   (`@allocated == 0`) and `test/test_type_stability.jl` (`@inferred` + JET
   `@test_opt`). Every phase runs the full suite; phases 5–7 additionally
   require the extended fixtures from phase 0.
4. **No new `@inbounds`.** Two existing ones are removed in §7.3 because they
   guard provably-literal indices the compiler already elides.

## Verification protocol

Run after every commit, not just at the end of a phase:

```julia
# in the MCP session (Revise picks up src/ edits automatically)
using Pkg; Pkg.test()
```

For the metaprogramming phases (5–6), the acceptance criterion is narrower and
must be checked *per converted function*, not just at suite level:

```julia
@allocated generate_equations(c) == 0          # for every fixture composite
@inferred  compute_residual(c, x, vars, others)
JET.@test_opt target_modules=(RheologyCalculator,) compute_residual(c, x, vars, others)
```

If a conversion regresses allocations or inference on *any* fixture, revert that
single conversion and leave the `@generated` in place with a comment stating the
constraint it satisfies. Partial completion of phase 6 is an acceptable outcome;
a suite that allocates is not.

## Status

- **Done.** Phase 0. `test_allocations.jl` now covers seven composites and
  `test_type_stability.jl` five, including the Kelvin-Voigt, generalized
  Maxwell, plastic, and volumetric topologies.
- **Done.** Phase 1.1 (`LTPViscosity` default argument), Phase 1.2A
  (delegating viscosity fallbacks), Phase 1.2B (fail-fast on a nonpositive KV
  aggregate), Phase 1.3 (`correct_xnorm` shape validation), and Phase 1.5
  (history-tuple bounds), all with regression coverage.
- **Done.** Phase 2. Every item deleted; `src/` is down from 4265 to 4011
  lines. §1.4 is settled with it — both functions that would have thrown are
  gone.
- **Done.** Phase 3, including a docs build that had been failing since the
  `NonConvergenceError` docstring was added.
- **Pending.** Phases 4–8.
- **Done.** `CompositeModel` removed (export, struct, and Documenter entry) — see §9.
- **Decided.** The tensor helpers become public API — see §4.4.

### Implementation log

- Added `test/test_simplification_plan.jl` covering the LTP default, delegated
  viscosity fallbacks, and valid/mismatched `xnorm0` shapes.
- Replaced the obfuscated `correct_xnorm` identity with matching-`SVector`
  dispatch plus explicit `DimensionMismatch`/`MethodError` behavior.
- Removed the obsolete `elastic_correction=false` keyword from
  `examples/elastic_corrections_playground.jl` as part of the preceding audit;
  the solver now applies the correction automatically.
- Added the four phase-0 fixtures to both `test_allocations.jl` and
  `test_type_stability.jl`. All four are allocation-free and JET-clean today, so
  they enter as a baseline rather than as a bug report.
- Added `_checked_η_KV` and routed the three `_η_KV` call sites through it
  (§1.2B). The `isfinite` clause the plan proposed is *not* included: a plastic
  leaf in parallel with a spring gives `η_KV = Inf` and a correctly zero
  correction, which is a legitimate topology. `η_KV > 0` alone rejects zero,
  negative, and `NaN`.
- Fixed §1.5, surfaced by the phase-0 Kelvin-Voigt fixture. Verified that the
  offending call returns a value with the `@inbounds` in place and raises a
  `BoundsError` without it, under default bounds checking.
- Verification after this batch: `Pkg.test()` passes (28 allocation tests, 35
  type-stability tests, 17 simplification-plan regressions), and
  `test/runtests.jl --allocations-only` passes without `--check-bounds=yes`.
- Phase 2 kept `global_series_functions` (called by `generate_equations`),
  `cpad` (called by `print_rheology_matrix`), and `compute_Q(::DruckerPrager)`
  (part of the plastic-model interface the examples use); the rest of §2.2–§2.7
  is deleted. `Pkg.test()` still passes. The `docs/make.jl` vitepress stage
  fails on this machine before and after the deletions, so it says nothing
  about them; Documenter itself reports no dangling `@docs` reference.

---

## Phase 0 — Widen the safety net (do this first)

*Correctness* coverage of the KV/generalized-Maxwell machinery is already good:
[test/test_strain_rate_correction.jl](test/test_strain_rate_correction.jl) is 329
lines with unit tests on `count_elastic`, `_n_elastic_in_parallel`, `_iselastic`,
`_η_eff_maxwell`, `_η_eff_elastic`, `_η_KV`, and `_weighted_backstress`
individually, plus eight end-to-end correction topologies. Phases 5–6 are
substantially de-risked by it, and it should be the first thing to run after any
change to `strain_rate_correction.jl`.

The gap is *performance* coverage. `test_allocations.jl` covers three composites,
all purely viscous/elastic and at most one nesting level.
`test_type_stability.jl` covers exactly one composite,
`SeriesModel(Elasticity, LinearViscosity)`. So a refactor that keeps every value
correct while quietly introducing an inference cliff or an allocation in the
KV/Maxwell or plastic paths passes CI today.

**Action.** Extend the `test_allocations.jl` and `test_type_stability.jl` fixture
sets before touching any source, adding:

- a Kelvin-Voigt branch: `SeriesModel(elastic, ParallelModel(elastic, viscous))`
- a generalized Maxwell branch:
  `SeriesModel(viscous, ParallelModel(SeriesModel(viscous, elastic), viscous))`
  — this is the only topology that exercises `_η_eff_maxwell`, `_η_eff_elastic`,
  and the `η_star = η_eff_M/η_el` weighting.
- a plastic composite (`DruckerPrager` in series) — exercises `compute_lambda`,
  `add_child`'s zero-returning specializations, and `subtract_parent`'s lambda
  exceptions.
- a volumetric composite (`Elasticity` + `BulkViscosity`) — exercises the
  `compute_pressure` / `compute_volumetric_strain_rate` equation elimination
  path.

**Why first:** phases 5–7 are mechanical, and `test_strain_rate_correction.jl`
will catch a wrong *number*. What it cannot catch is the characteristic failure
of these refactors — code that stays correct but stops inferring, because the
compiler no longer unrolls what a `@generated` used to unroll by construction.
Only an `@allocated`/`@inferred` fixture on the affected topology catches that.

**Risk:** none (test-only).
**Estimated size:** +60 lines of test code.

---

## Phase 1 — Fix the latent defects

Four defects surfaced during the read. Three are confirmed by reproduction on
this machine; the fourth is a robustness gap of the same family. They are fixed
first because they are independent of every refactor below, and because §1.2 in
particular is a *silent wrong-answer* path — the sharpest kind of fail-fast
violation.

Each fix lands with a regression test in the same commit.

### 1.1 `LTPViscosity` throws `UndefVarError` on its own default argument

[LTP_viscosity.jl:31](src/rheology/viscous/nonNewtonian/LTP_viscosity.jl#L31):

```julia
@inline compute_viscosity_parallel(r::LTPViscosity; τ = 00e, kwargs...) = ...
```

Julia parses `00e` as `0 * e`. There is no `e` in scope, so the default is a
reference to an undefined variable.

**Reproduction** (verified):

```julia
julia> compute_viscosity_parallel(LTPViscosity(6.2e-13, 76.0, 1.8e9, 3.4e9))
ERROR: UndefVarError: `e` not defined in `RheologyCalculator.RheologyModels`
```

Every call that does not pass an explicit `τ` throws. It has gone unnoticed
because `LTPViscosity` defines no `parallel_state_functions`, and the one
internal caller — `compute_viscosity_parallel(::NTuple, ε, τ, others)` in
`strain_rate_correction.jl` — always merges a `τ` into its argument tuple.

**Fix.** `τ = 0e0`, matching the two lines above it. Note this method is *not*
redundant and so is not deleted by §4.3: `LTPViscosity`'s parallel effective
viscosity `τ / 2ε̇(τ)` genuinely differs from its series form.

**Regression test.** Call all three `compute_viscosity*` methods of every
exported element with no keyword arguments and assert none throws. This is three
lines and would have caught the defect at introduction.

### 1.2 A missing viscosity method silently halves the effective viscosity

`src/core/state_functions.jl` gives every state function, including
`compute_viscosity`, `compute_viscosity_series`, and `compute_viscosity_parallel`,
a blanket `AbstractRheology` fallback returning `0.0`. For the *residual* state
functions that is correct — an element that does not participate in an equation
contributes nothing to the sum. For the *viscosity* aggregates it is not: `_η_KV`
sums them, so a zero contribution is indistinguishable from "this element has no
stiffness."

**Reproduction** (verified). An element that defines `compute_viscosity` but not
`compute_viscosity_series` — exactly what happens when a new element file omits
an import, see §4.2:

```
η_KV, all elements complete : 2.0e20
η_KV, one element forgetful : 1.0e20     <- silently halved, no error
```

And an element that defines none of the three:

```
η_KV, branch with no viscosity methods at all: 0.0
  -> ws / (2*η_KV) with ws=1e6 gives: Inf   (silent, no error)
```

The first case is the dangerous one: the KV correction comes back a factor of two
wrong, the Newton iteration converges happily, and the result is simply incorrect.

**Fix, part A — make the fallback delegate rather than zero.** In
`src/core/state_functions.jl`, exclude the two composition-specific viscosities
from the blanket `0.0` loop and define them as delegating to `compute_viscosity`:

```julia
# Composition-specific effective viscosities default to the element's own
# effective viscosity. Unlike the residual state functions, these are summed by
# _η_KV, where a zero contribution is indistinguishable from a real one — so the
# fallback must delegate rather than return zero.
@inline compute_viscosity_series(r::AbstractRheology; kwargs...)   = compute_viscosity(r; kwargs...)
@inline compute_viscosity_parallel(r::AbstractRheology; kwargs...) = compute_viscosity(r; kwargs...)
```

This resolves the first case exactly, and it is behavior-preserving for every
element currently in the package: all six that define `compute_viscosity_series`
define it identically to `compute_viscosity` (verified textually), and elements
that define none still reach `compute_viscosity`'s own `0.0` fallback.

Part A also removes ~12 lines of boilerplate, which is why §4.3 is now only the
deletion half of this change.

**Fix, part B — fail fast on a zero aggregate.** Part A cannot rescue an element
that defines no viscosity at all, nor a caller that omits `dt` (which makes an
elastic leaf's `G*dt` zero). Both currently produce `Inf` corrections. `_η_KV` is
only ever called on a branch that has already been established to contain an
elastic element, so a zero or non-finite aggregate is unambiguously a bug, never
a legitimate state. Add the check at the point of use in
`_kv_branch_correction`, `_branch_elastic_stress`, and
`_kv_implicit_branch_correction_scalar`:

```julia
isfinite(η_KV) && η_KV > 0 || throw(ArgumentError(
    "effective Kelvin-Voigt viscosity of this parallel branch is $η_KV; every " *
    "element in a branch carrying elastic history must define compute_viscosity, " *
    "and `others` must supply a nonzero `dt`"))
```

A plastic-only branch gives `η_KV = Inf` (`compute_viscosity(::DruckerPrager) = Inf`),
which the `isfinite` clause would reject — confirm against the phase-0 plastic
fixture whether such a branch can reach `_η_KV` at all, and if so relax to
`η_KV > 0` only.

### 1.3 `correct_xnorm` accepts anything and validates nothing

[normalize_x.jl:40](src/equation_system/normalize_x.jl#L40):

```julia
@inline correct_xnorm(::SVector{N,T}, xnorm) where {N, T} = (T; xnorm)
```

`(T; xnorm)` is a block that evaluates the type parameter `T`, discards it, and
returns `xnorm` — an obfuscated identity function. It performs no check that
`xnorm` is a vector, or that its length matches `x`.

**Reproduction** (verified):

```julia
julia> correct_xnorm(SA[1.0, 2.0], "garbage")
"garbage"
```

A mismatched or malformed `xnorm0` passed to `solve` therefore fails later inside
`mynorm`, with an error that points at the norm rather than at the caller's
argument.

**Fix.** Constrain the signature and check the length:

```julia
@inline correct_xnorm(::SVector{N}, xnorm::SVector{N}) where {N} = xnorm
@inline correct_xnorm(::SVector{N}, ::SVector{M}) where {N, M} = throw(DimensionMismatch(
    "xnorm0 has $M entries but the solver vector has $N"))
@inline correct_xnorm(::SVector{N, T}, ::Nothing) where {N, T} = @SVector ones(T, N)
```

Every current caller passes either `nothing` or the output of `normalisation_x`,
which is an `SVector` of matching length, so no working call changes.

### 1.4 Two dead functions would throw if they were ever called

Not fixes so much as evidence for the deletions in §2. Recording them here so the
reasoning is not lost:

- `global_parallel_functions` and `local_parallel_functions`
  ([composite.jl](src/core/composite.jl)) both call `count_parallel_elements`,
  **which is not defined anywhere in the package**. Either would throw
  `UndefVarError` on first call.
- `create_rheology_string(::ParallelModel)`
  ([print_rheology.jl:195](src/display/print_rheology.jl#L195)) reads
  `rheo_Parallel.elements`; `ParallelModel` has no such field.

Both are deleted in §2 rather than repaired — nothing calls them, and neither has
a specification to repair them against.

**Risk:** §1.1 and §1.3 are contained. §1.2 part A is behavior-preserving for
every element in the package but changes the contract for *downstream* elements
defined outside it — an external element relying on `compute_viscosity_series`
returning `0.0` while `compute_viscosity` returns nonzero would change behavior.
That combination is already a bug, but mention it in the release notes.
**Reduction:** ~12 lines (from §1.2 part A).

### 1.5 A short history tuple reads arbitrary memory

[equations.jl:356](src/equation_system/equations.jl#L356) wrapped the whole
per-field extraction in `@inbounds`:

```julia
Base.@ntuple $N i -> @inbounds _extract_local_kwargs(vals_args[i], keys_args[i], keys_hist, n)
```

`_extract_local_kwargs` indexes a history tuple by the element number,
`vals_args[n]`. `history_kwargs(::AbstractElasticity) = (:τ0, :P0)`, so every
elastic element claims *both* `τ0` and `P0` — including `IncompressibleElasticity`,
which never reads `P0`. A composite with two springs therefore requires two
entries in each of `τ0` and `P0`, and an `others` carrying `P0 = (0.0,)` reads
`P0[2]` out of range. Under `@inbounds` that is silent undefined behavior; the
value substituted for the missing entry is whatever memory follows the tuple.

**Reproduction** (verified): `SeriesModel(elastic, ParallelModel(elastic, viscous))`
with `others = (; dt, τ0 = (0.0, 0.0), P0 = (0.0,))` returns an initial guess in
an ordinary session and throws `BoundsError` under `--check-bounds=yes` — which
is what `Pkg.test` runs, and how the phase-0 Kelvin-Voigt fixture surfaced it.

**Fix.** Drop the `@inbounds`, and with it the now-pointless
`Base.@propagate_inbounds` on `_extract_local_kwargs`. Julia's own bounds check
is the whole fix: the same call that returned a value now raises a `BoundsError`
under default settings, and no bespoke guard is needed to say so.

This is the `@inbounds` removal §7.3 already scheduled for
[equations.jl:356](src/equation_system/equations.jl#L356); it lands here instead
because it is a silent-wrong-answer path, not a cosmetic cleanup. The five
`@inbounds` in `composite.jl` remain scheduled for §5.2.

**Risk:** none for correct callers — the check only fires where the old code
read out of range. Allocation-free (`Pkg.test` and the non-`check-bounds`
allocation run both report zero).

---

## Phase 2 — Delete unreachable code

Roughly 250 lines of `src/` (~6%) is never called from `src/`, `ext/`, `test/`,
`examples/`, `docs/`, `prototypes/`, or `2D/`. Verified by name-grep across the
whole repository.

### 2.1 Two dead `add_global_equations` methods

`src/equation_system/equations.jl:245-267`. Both `@generated` methods take 9
positional arguments. The only call site,
[equations.jl:80](src/equation_system/equations.jl#L80), passes 11 and is
matched by the plain method at
[equations.jl:269](src/equation_system/equations.jl#L269). The two `@generated`
methods are unreachable.

Delete both (23 lines). This also removes 2 of the 63 `@generated` functions
before phase 6 has to consider them.

### 2.2 Dead composite helpers

`src/core/composite.jl`:

- `global_parallel_functions` — no callers, and calls the undefined
  `count_parallel_elements` (§1.4).
- `local_parallel_functions` — same.
- `local_series_functions` — no callers.
- `hasbranches` — no callers.

### 2.3 Dead argument helpers

`src/core/kwargs.jl` — `augment_args`, `update_args2`, `residual_kwargs`,
`all_differentiable_kwargs`, `split_args` (~74 lines, all with docstrings). None
has a caller. `residual_kwargs` is referenced from
[docs/src/api.md](docs/src/api.md); remove that entry in the same commit.

This leaves `kwargs.jl` at `history_kwargs` + `differentiable_kwargs` — the two
that equation generation actually uses — which makes the file's purpose legible
for the first time.

### 2.4 Dead rheology-type helpers

`src/core/rheology_types.jl` — `length_state_functions` and
`get_unique_state_functions`. ~31 lines. `flatten_repeated_functions` is the live
deduplication path; `get_unique_state_functions` is a `Symbol`-dispatched wrapper
around it that nothing calls.

### 2.5 Dead display code

`src/display/print_rheology.jl`:

- `create_rheology_string` (22 lines) — no callers, and broken (§1.4).
- `create_parallel_str` (42 lines) — no callers.
- `InverseCreepLaw` (11 lines) — an `AbstractRheology` subtype with no state
  functions and no constructor call anywhere; it is also unrelated to printing
  and simply lives in the wrong file.

### 2.6 Dead predicate

`is_eq_elastic` ([post_calculations.jl:4-5](src/post_processing/post_calculations.jl#L4-L5)) — no callers.

### 2.7 Commented-out code

`src/core/others.jl:47-59` (the `harmonic_average` block) and the commented
alternative branches in `Drucker_Prager.jl` (lines 38-42, 46-50, 68, 87, 92).
Git history preserves these; the source should not.

Where a commented line records an unresolved design question rather than a
former implementation — e.g. `Drucker_Prager.jl:37` "we need to check whether
this allocates" — convert it to a GitHub issue and reference the issue number,
or delete it.

**Risk:** minimal — everything here is provably unreachable.
**Reduction:** 254 lines. **Done.**

---

## Phase 3 — Fix stale and incorrect documentation

The comments in `src/` are unusually dense, which is good, but a significant
fraction point at files that do not exist. A reader who follows them wastes time
and then stops trusting the rest.

### 3.1 Dangling file references (14 occurrences)

| Reference | Occurrences | Status |
|---|---|---|
| `tensor_reduction.typ` | 10 | **No such file in the repository** |
| `RheologyDefinitions.jl` | 2 | **No such file** |
| `CLAUDE.md` | 1 | **No such file** |
| `src/strain_rate_correction.jl` | 3 | Wrong path — the file is `src/post_processing/strain_rate_correction.jl` |

`tensor_reduction.typ` is the worst case: ten comments defer the *entire
derivation* of the KV/generalized-Maxwell correction to a document that is not
in the repository. The math it refers to is the least obvious code in the
package.

**Action:**

- Fix the three `src/strain_rate_correction.jl` paths (also in
  [docs/src/strain_rate_correction.md:293](docs/src/strain_rate_correction.md#L293)).
- Delete the `RheologyDefinitions.jl` and `CLAUDE.md` references; state the
  invariant directly instead. For `equations.jl:102`, the substance is already
  in the comment ("equation position is assumed to equal `.self`") — just drop
  the citation.
- For `tensor_reduction.typ`, choose one: **(a)** commit the document to the
  repository (e.g. `docs/src/derivations/`) and keep the references, or **(b)**
  move its conclusions — equations (\*) for `ε_eff` and the `η_star` definitions,
  which the comments already restate — into
  [docs/src/strain_rate_correction.md](docs/src/strain_rate_correction.md) and
  point there. Option (a) is better if the derivation exists somewhere; the
  comments are currently the only record of where the formulas come from.

### 3.2 Wrong provenance comment

[mod_Cam_Clay.jl:1-3](src/rheology/plastic/critical_state/mod_Cam_Clay.jl#L1-L3)
carries the header "This implements the mode1/mode2 plasticity model proposed in
Popov et al. (2025)" — copy-pasted from `Drucker_Prager_cap.jl`. The docstring
two lines below correctly cites de Souza Neto. Delete the copied header.

### 3.3 Misleading physics comment

[strain_rate_correction.jl:137-138](src/post_processing/strain_rate_correction.jl#L137-L138)
labels `compute_viscosity_parallel` over a tuple as "Harmonic-mean effective
viscosity of a tuple of elements in *parallel* (stresses add → strain rates must
be consistent → harmonic mean)". Elements in parallel add stresses at common
strain rate, which is an *arithmetic* sum of viscosities. The code computes a
harmonic mean in both the series and parallel methods.

**Do not change the code.** Establish first whether the harmonic form is
intentional (these are aggregation helpers for initial guesses and corrections,
not constitutive law) and then write a comment that states what is actually
true. If it is a genuine bug, that is a separate PR with its own regression
test — unlike the phase-1 defects, this one has no reproduction yet.

### 3.4 Duplicate Documenter entry

`RheologyCalculator.compute_residual` appears twice in
[docs/src/api.md](docs/src/api.md) — once under "Solver", once under
"Internals". Documenter normally errors on a duplicated docstring; this is
silently swallowed by `warnonly` (§8.2). Remove one.

### 3.5 Shadowed argument in `ModCamClay`

[mod_Cam_Clay.jl:59-68](src/rheology/plastic/critical_state/mod_Cam_Clay.jl#L59-L68)
and the matching `compute_Q`:

```julia
function compute_F(r::ModCamClay, τII, P)
    (; M, r, β, Pt) = r    # `r` the rheology is rebound to `r` the radius
```

This works — the right-hand side is evaluated before the destructuring binds —
but every subsequent `r` in the function means the radius, not the element. Rename
the argument to `p` (or the field access to `p.r`) while §5.1 is touching these
files anyway.

**Risk:** none (comments, docs, and a local rename). **Done.**

### 3.6 Findings from carrying out this phase

- `docs/derivations/tensor_reduction.typ` is the restored derivation: it was
  removed in d7613d3 and recovered from 07a9513. `GenMaxwell.typ` and
  `eff_strain.jl` went in the same commit and are still only in history; no
  comment cites them.
- §3.4 was not what broke the documentation build. A duplicated `@docs` entry
  builds fine; what failed the vitepress render was `NonConvergenceError`
  having a docstring in the module and no entry in any `@docs` block. Verified
  by toggling each change independently. The build has been failing since that
  docstring was introduced, which is the concrete cost of the `warnonly`
  setting §8.2 addresses.
- §3.3 needed no comment at all. The two four-argument tuple reductions it
  describes have no caller at that arity anywhere in the repository, so they
  are deleted rather than re-documented, and the harmonic-versus-arithmetic
  question does not arise. `_η_KV` (arithmetic, parallel) and
  `_η_eff_maxwell` (harmonic, series) are the live paths and are each correct
  for their composition.

---

## Phase 4 — Collapse per-element boilerplate

This is the highest ratio of lines-removed to risk-taken in the plan.

### 4.1 The cap-family plastic models — ~135 lines

`DruckerPragerCap`, `ModCamClay`, `Hyperbolic`, and `Golchin` each define the
same eight methods. Verified by textual diff after normalizing the type name:
**the method blocks are character-for-character identical across all four
files.** They differ only in the struct, its constructor, `compute_F`, and
`compute_Q`.

The shared block, in every one of the four files:

```julia
@inline series_state_functions(::M)   = (compute_strain_rate, compute_lambda, compute_volumetric_strain_rate)
@inline parallel_state_functions(::M) = compute_stress, compute_pressure, compute_lambda,
                                        compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate
@inline _isvolumetric(::M) = true
compute_strain_rate, compute_volumetric_strain_rate, compute_lambda,
compute_stress, compute_pressure,
compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate
```

**How.** Introduce an abstract supertype in `src/core/rheology_types.jl`:

```julia
"""
    AbstractCapPlasticity <: AbstractPlasticity

Plastic elements whose yield surface is closed in pressure and whose flow rule
is derived from a potential by automatic differentiation. Subtypes provide
`compute_F(r, τII, P)` and `compute_Q(r, τII, P)`; the state functions, the
Duvaut-Lions regularization, and the flow rule are inherited.
"""
abstract type AbstractCapPlasticity <: AbstractPlasticity end
```

Write the eight shared methods **once**, against `::AbstractCapPlasticity`, in a
new `src/rheology/plastic/cap_plasticity.jl`. Each of the four model files then
retains only: its citation header, its struct, its keyword constructor, and its
`compute_F`/`compute_Q`. Each drops from ~100 to ~55 lines, and adding a fifth
cap model becomes a ~40-line file instead of a ~100-line copy-paste.

`DruckerPrager` (the plain one) is **not** part of this family and must stay
separate: its `compute_lambda` uses a unit Lagrange multiplier rather than
`η_vp`, its `compute_plastic_strain_rate` is hardcoded to `λ - ε` rather than
AD-differentiated, and `_isvolumetric` is `false`.

**Verification.** Dispatch on an abstract supertype is the same mechanism the
package already uses for `AbstractElasticity`/`AbstractViscosity` in
`history_kwargs`, so type stability is unaffected. Confirm with the existing
`test_ModCamClay.jl`, `test_VEPCap.jl`, `test_Hyperbolic.jl` and the phase-0
plastic fixture.

**Risk:** low. **Reduction:** ~135 lines, and each future model ~60 lines cheaper.

### 4.2 The repeated import headers — ~39 lines

Twelve of the sixteen files under `src/rheology/` open with 2–5 lines of
`import ..RheologyCalculator: ...` (45 such lines total), listing overlapping
subsets of the same ~15 names.

Beyond the line count this is the mechanism behind §1.2: a new element file that
*forgets* to import, say, `compute_viscosity_series` silently defines a **new
local function** of that name in `RheologyModels` rather than extending the core
one. The core generic then keeps its fallback for that element. Phase 1 makes the
consequence benign (the fallback now delegates to `compute_viscosity`); this
phase removes the opportunity.

**How.** Move a single import block to `src/RheologyModels.jl`, immediately
before the `include`s:

```julia
# Names extended by the element files below. Importing them here rather than
# per-file ensures every `include`d method extends the core generic instead of
# defining a same-named function local to this module.
import ..RheologyCalculator: series_state_functions, parallel_state_functions, _isvolumetric,
    compute_strain_rate, compute_stress, compute_pressure, compute_volumetric_strain_rate,
    compute_plastic_strain_rate, compute_plastic_stress, compute_volumetric_plastic_strain_rate,
    compute_lambda, compute_lambda_parallel,
    compute_viscosity, compute_viscosity_series, compute_viscosity_parallel,
    isvolumetric
```

Delete all 45 per-file import lines. The `import ForwardDiff: ForwardDiff` lines
in the plastic files can also move up, since `RheologyModels.jl` already imports
it at line 11.

**Verification.** Add `ExplicitImports.jl` to the test suite (the repo has a
`freshen-explicit-imports` workflow available) so this class of mistake is caught
mechanically rather than by review. Run `Pkg.test()`, and diff
`methods(compute_viscosity_series)` before and after — the method count must be
unchanged.

**Risk:** low, but this is the one refactor where a mistake is *quiet*. Own commit.

### 4.3 Delete the redundant `compute_viscosity` methods — ~12 lines

Phase 1.2 part A installs the delegating fallbacks. This item is the deletion
half: remove every `compute_viscosity_series` definition (all six are textually
identical to the element's `compute_viscosity`) and the three
`compute_viscosity_parallel` definitions that merely restate `compute_viscosity`
(`LinearViscosity`, `Elasticity`, `IncompressibleElasticity`, `DruckerPrager`).

Keep the three genuinely different `_parallel` methods on the nonlinear laws
(`PowerLawViscosity`, `DiffusionCreep`, `LTPViscosity`) — these compute
`τ / 2ε̇(τ)` rather than `σ(ε̇) / 2ε̇`, which is a different quantity.

### 4.4 The 2D/3D tensor helpers — ~30 lines, and a public API

[src/utils/tensor_helpers.jl](src/utils/tensor_helpers.jl) is five pairs of
functions that differ only in tuple length: `second_invariant_{2D,3D}`,
`tensor_strain_rate_{2D,3D}`, `vars_{2D,3D}`, `zero_stress_tensor_{2D,3D}`,
`stress_tensor_from_invariant_{2D,3D}`, `elastic_stress_history_{2D,3D}`.
`second_invariant_2D`/`_3D` additionally duplicate `second_invariant` from
`strain_rate_correction.jl`.

**(a) Internal deduplication.** Make the bodies dimension-generic by dispatching
on `NTuple{3}` / `NTuple{6}` and delegating to the existing `second_invariant`.
The `_2D`/`_3D` names are used in **31 test and example files** and must be kept
as thin wrappers — this is an internal cleanup, not a rename.

**(b) Make them public — decided.** These helpers are used by essentially every
test and every example, yet they are currently neither exported nor documented,
so `test/runtests.jl` reaches them with an explicit
`import RheologyCalculator.RheologyModels: second_invariant_2D, ...` and a
comment apologizing for it. Every user who copies an example hits the same wall.

Action:

1. Add to the `export` list in `src/RheologyModels.jl`:
   `second_invariant_2D`, `second_invariant_3D`, `tensor_strain_rate_2D`,
   `tensor_strain_rate_3D`, `vars_2D`, `vars_3D`, `zero_stress_tensor_2D`,
   `zero_stress_tensor_3D`, `stress_tensor_from_invariant_2D`,
   `stress_tensor_from_invariant_3D`, `elastic_stress_history_2D`,
   `elastic_stress_history_3D`.
2. Give each a docstring — they currently have none. State the Voigt component
   order explicitly (`(xx, yy, xy)` and `(xx, yy, zz, yz, xz, xy)`), which is
   currently only inferable from the arithmetic in `second_invariant`.
3. Add a "Tensor Helpers" section to [docs/src/api.md](docs/src/api.md).
4. Delete the `import RheologyCalculator.RheologyModels: ...` line and its
   apologetic comment from `test/runtests.jl`.
5. Update the `RheologyModels.jl` header comment, which currently says advanced
   models are "intentionally NOT exported" — clarify that this refers to the
   material models, not the tensor utilities.

Since this widens the exported surface, `Aqua.test_undefined_exports` (already in
the suite) covers it, and the new docstrings make `checkdocs` meaningful for them.

**Risk:** (a) low; (b) additive — no existing name changes meaning.

---

## Phase 5 — Deduplicate the structural helpers

These are pairs and triples of functions with identical shape, differing only in
which of `series`/`parallel` or `tensor`/`scalar` they operate on. Each is a
place where a bug fix has to be applied twice and might not be.

### 5.1 `series_*` / `parallel_*` leaf and branch extraction

[composite.jl](src/core/composite.jl) — four five-method blocks (`series_leafs`,
`parallel_leafs`, `series_branches`, `parallel_branches`) with identical
recursion structure; only the "which composite type is a branch" predicate
differs (`ParallelModel` vs `SeriesModel`).

**How.** Parametrize on the branch type and generate the four with a small
`for`-loop over `(:series => ParallelModel, :parallel => SeriesModel)`, in the
style the file already uses for `_local_series_state_functions` and
`_local_parallel_state_functions`. 23 lines become ~10.

### 5.2 `local_`/`global_` × `series_`/`parallel_` state functions

Four `@generated` tuple-filters in `composite.jl` with identical bodies. Same
treatment: one shared implementation, four thin dispatching wrappers. Removes 4
of the 63 `@generated`.

### 5.3 `_estimate_initial_value_harm` / `_arith`

[initial_guess.jl:168-224](src/equation_system/initial_guess.jl#L168) — two
`@generated` functions with ~40 lines of docstring each and bodies that differ
in exactly one line (`sum_vals += safe_inv(val)` vs `sum_vals += val`) plus the
return (`safe_inv_one(sum_vals)` vs `sum_vals`).

**How.** One implementation taking the accumulation and finalization as
arguments; the two entry points become one-liners. The two long docstrings
collapse to one that documents both estimates in a two-row table.

### 5.4 `_weighted_backstress` / `_weighted_backstress_scalar`

[strain_rate_correction.jl:484](src/post_processing/strain_rate_correction.jl#L484)
and [660](src/post_processing/strain_rate_correction.jl#L660). Structurally
identical `@generated` functions — the comment on the second already says so
("Structurally identical to `_weighted_backstress` but accumulates a scalar").
They differ only in whether the accumulator starts at `ε .* 0` or `0.0` and
whether each `τ0` entry passes through `second_invariant_value`.

**How.** One `@generated` parametrized by an accumulator seed and a per-entry
reduction; two one-line entry points. This is worth doing specifically because
these two must stay in agreement — the implicit residual correction and the
pre-solve tensor correction have to compute the same weighting, and right now
nothing enforces that beyond both files being edited together.

### 5.5 The three per-branch offset walkers

`_kv_corrections`, `_kv_corrections_elastic_stress`, and
`_kv_implicit_corrections_scalar` all begin with the same three lines:

```julia
foreach(_assert_kv_nesting_supported, branches.parameters)
counts  = [_n_elastic_in_parallel(branches.parameters[i]) for i in 1:Nb]
offsets = cumsum([0; counts[1:(end - 1)]])
```

**How.** Extract a compile-time helper `_branch_tau0_offsets(branches_type)`
returning the validated offset vector, called by all three. The offset
convention (τ0 ordered leafs-first, then branches, matching
`global_eltype_numbering`) is then stated in exactly one place instead of
restated in three comment blocks.

Similarly, the `findall(Ti -> Ti <: AbstractElasticity, ...)` + `sub_elastic_pos`
scan appears **five times** (`_weighted_backstress`,
`_weighted_backstress_scalar`, `_branch_elastic_info`,
`_direct_leaf_correction_scalar`, `_n_elastic_in_parallel`). Extract one
`_elastic_source_positions(leafs_type, subs_type)`.

### 5.6 `_iselastic` over tuples

Two identical `@generated` functions, one for `NTuple{N,AbstractRheology}` and
one for `NTuple{N,AbstractCompositeModel}`, plus a `Tuple{}` tiebreaker whose
comment explains that the split exists only to resolve the `N=0` ambiguity. A
single method on `NTuple{N,Union{AbstractRheology,AbstractCompositeModel}}` plus
the `Tuple{}` method removes the duplicate and the tiebreaker comment with it.

### 5.7 `compute_viscosity_series` / `_parallel` over tuples

[strain_rate_correction.jl:139-166](src/post_processing/strain_rate_correction.jl#L139) —
two `@generated` functions with byte-identical bodies (see §3.3 on the comment).
Generate both from a `for` loop over the two names.

**Risk:** moderate. All of these are on inference-critical paths. Convert one at
a time; run the per-function `@inferred`/`@allocated` checks after each.
**Reduction:** ~120 lines, and 6–8 fewer `@generated` functions.

---

## Phase 6 — Centralize the metaprogramming

This is the single largest readability lever and also the one to approach most
carefully.

**The finding.** `src/` contains **63 `@generated` functions in 4265 lines** —
one every 68 lines. Classified by body:

| Pattern | Count | Verdict |
|---|---:|---|
| Plain tuple map: `Base.@ntuple $N i -> f(x[i])` | 40 | Replaceable by one shared helper |
| Accumulator: `Base.@nexprs $N i -> acc += …` | 14 | Replaceable by one shared helper |
| Type introspection (`.parameters`, `findall`, `cumsum`) | 9 | **Genuinely requires `@generated`** — keep |

Two of the 40 are the dead methods deleted in §2.1, and 4–8 more disappear in
phase 5, so the live target is roughly 30 maps and 12 accumulators.

**The recommendation: do not remove `@generated` — centralize it.**

The naive conversion — `Base.@ntuple $N i -> f(x[i])` → `ntuple(i -> f(x[i]), Val(N))`
— is tempting but wrong here. `Base.ntuple(f, ::Val{N})` is explicitly unrolled
only for `N ≤ 10`; beyond that it falls back to a path that does not infer for
heterogeneous results. Many of these tuples *are* heterogeneous — `eqs` is an
`NTuple{N,CompositeEquation}` whose elements have different type parameters, and
the state-function tuples are tuples of distinct function types. A composite with
more than ten equations is entirely plausible (a VEVP model with several parallel
branches), and the failure mode would be a silent inference cliff that only shows
up as allocations on large models.

Instead, define **two** helpers, in a new `src/core/tuple_utils.jl`, each written
as recursion over `first`/`Base.tail` (which infers and unrolls for arbitrary `N`
without any `@generated` at all):

```julia
"""
    maptuple(f, t::Tuple)

Apply `f` to each element of `t` and return the results as a tuple. Unrolls and
infers for tuples of mixed element types and arbitrary length, which
`ntuple(f, Val(N))` does not do beyond `N = 10`.
"""
@inline maptuple(f::F, t::Tuple) where {F} = (f(first(t)), maptuple(f, Base.tail(t))...)
@inline maptuple(::F, ::Tuple{}) where {F} = ()

"""
    foldtuple(op, init, f, t::Tuple)

Left-fold `op` over `f` applied to each element of `t`, starting from `init`.
"""
@inline foldtuple(op::O, init, f::F, t::Tuple) where {O, F} =
    foldtuple(op, op(init, f(first(t))), f, Base.tail(t))
@inline foldtuple(::O, acc, ::F, ::Tuple{}) where {O, F} = acc
```

Then rewrite the ~30 map sites and ~12 accumulator sites as calls. A
representative before/after — `evaluate_state_functions`
([equations.jl:374](src/equation_system/equations.jl#L374)):

```julia
# before: 6 lines, a @generated function, a quote block, an @inline inside the quote
@generated function evaluate_state_functions(eqs::NTuple{N, CompositeEquation}, args, others) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> evaluate_state_function(eqs[i], args[i], others)
    end
end

# after: 1 line
@inline evaluate_state_functions(eqs::Tuple, args, others) =
    maptuple2((eq, a) -> evaluate_state_function(eq, a, others), eqs, args)
```

Several sites zip `eqs` against `args` or `x` rather than mapping one tuple, so
add a third helper `maptuple2(f, t1, t2)` walking two tuples in lockstep — same
recursion, and it avoids index arithmetic entirely.

**What this buys.** The metaprogramming that a reader must understand drops from
63 sites to 3 helpers plus 9 genuinely-type-introspecting functions. The 9 that
remain are the ones where `@generated` is actually load-bearing — they inspect
`.parameters` to find elastic elements and bake τ0 indices in as literals — and
with 54 neighbors gone, they stand out as the deliberate choices they are, each
already carrying a good explanatory comment.

**How to sequence it.** One function per commit, in this order (easiest and
best-covered first):

1. `src/core/rheology_types.jl` (3 sites) and `src/core/others.jl` (2) — small,
   directly covered by `test_equation_graphs.jl`.
2. `src/core/composite.jl` (5, minus those already merged in §5.2).
3. `src/equation_system/initial_guess.jl` (5) and `normalize_x.jl` (1).
4. `src/equation_system/solver.jl` (3) — note `mynorm` and `max_feasible_step`
   walk `SVector`s, not tuples; they need the `SVector` form of the fold or can
   simply be left alone (they are 10 lines each and already clear).
5. `src/equation_system/equations.jl` (~18 after §2.1) — the hot path; run the
   full allocation and JET suite after **every single** conversion here.
6. `src/post_processing/post_calculations.jl` (4).

Stop wherever the numbers stop cooperating. The value is roughly linear in the
number of sites converted, so a partial phase 6 is still a win.

**Risk:** moderate-to-high, fully mitigated by the phase-0 fixtures and the
per-function acceptance criteria. This phase is the reason phase 0 exists.
**Reduction:** ~150 lines, and a step change in how much Julia arcana a
contributor needs to read the equation system.

---

## Phase 7 — `print_rheology.jl`

After §2.5 removes the dead half of the file, ~180 lines remain, with three
problems.

### 7.1 Three near-identical matrix builders

`print_rheology_matrix` has separate methods for `Tuple`, `ParallelModel`, and
`SeriesModel`. All three run the same "place sub-blocks into a matrix, fill
unassigned entries with empty strings, trim" loop; they differ in whether a
sub-block advances the row or the column cursor, and in whether the result gets
bracketed.

**How.** One builder taking the growth direction (`:row` for parallel, `:column`
for series) and a post-processing step. The `Tuple` and `SeriesModel` methods
are already nearly identical to each other — those two collapse with no
behavioral question at all.

### 7.2 Fixed 40×40 buffers silently truncate

All three methods allocate `Matrix{String}(undef, 40, 40)` and then write into
it with no bounds discipline. A composite with more than 40 elements in a row
writes out of range — `A[i:(i+si[1]-1), j]` on an oversized block throws
`BoundsError` at best, and the trimming step `A = A[1:i_max, 1:j]` silently drops
content in other paths.

Given a stated preference for failing fast: size the matrix from the actual
element count (available from `global_eltype_numbering`), or keep the buffer and
add an explicit check with a clear message. Either is better than 40.

Related: `length_str_no_colors(str) = textwidth(...) + 1 * 6 + 2` — `1 * 6 + 2` is
an unexplained constant, and `1 *` is inert. Name it or explain it.

### 7.3 Minor

- `import Base.show` at line 2 is redundant with the `function Base.show`
  definitions below. Drop the import.
- `Base.show(::IO, ::SeriesModel)` and `Base.show(::IO, ::ParallelModel)` have
  identical bodies — one method on `Union{SeriesModel, ParallelModel}`.
- The five `@inbounds` in `composite.jl` wrap `funs[i]` where `i` is a literal
  index generated by `Base.@ntuple $N` and `N == length(funs)` by the signature.
  The compiler elides these checks on its own; the annotation buys nothing and
  costs the guarantee that an indexing mistake in a future edit throws instead of
  reading arbitrary memory. Remove them as part of §5.2, which rewrites those
  functions anyway. The one in
  [equations.jl:356](src/equation_system/equations.jl#L356) is already gone —
  see §1.5.

**Risk:** low — this code has no numerical role. It is also currently untested;
add a smoke test that `repr(c)` runs without error for each phase-0 fixture,
which is cheap and would have caught the `.elements` bug in §1.4.

---

## Phase 8 — Tooling that keeps this from coming back

### 8.1 Two unused dependencies

`LaTeXStrings` and `TimerOutputs` are in `[deps]` and `[compat]` but appear
**zero times** in `src/`, `ext/`, and `test/`. Remove both from `Project.toml`,
run `Pkg.resolve()`.

They survived because `test/test_Aqua.jl` runs `test_deps_compat` but **not**
`Aqua.test_stale_deps`. Add it. Also note the entire Aqua testset is gated behind
`if VERSION ≤ v"1.12.3"`, so on the current 1.12.7 **none of it runs** — including
the ambiguity and piracy checks. Investigate why that gate exists (presumably a
transient Aqua/Julia incompatibility) and either narrow it to the specific
failing check or remove it.

### 8.2 Documenter swallows its own errors

[docs/make.jl:20](docs/make.jl#L20) sets `warnonly = Documenter.except(:footnote)`,
which downgrades *every* Documenter error except footnotes to a warning —
missing docstrings, broken cross-references, and the duplicate `compute_residual`
entry from §3.4 all pass silently. Combined with `checkdocs = :exports`, which
does not check the many documented-but-unexported internals, the docs build
cannot currently fail.

Tighten incrementally: fix the errors that surface, then narrow `warnonly` to the
specific categories that genuinely cannot be satisfied. This is the mechanism
that would have prevented §3.1's fourteen dangling references from accumulating.

### 8.3 Formatting

There is no `.JuliaFormatter.toml` and no format check in CI, and the style is
correspondingly inconsistent — `NTuple{N,Any}` vs `NTuple{N, Any}`, `where N` vs
`where {N}`, `1e-8` vs `1.0e-8` — sometimes within one file. Adopt `runic` (there
is a `freshen-runic` workflow available), commit the reformat as a standalone
commit, and add it to `.git-blame-ignore-revs`.

**Do this last**, after all the structural phases, so the reformat commit does not
collide with in-flight refactors.

---

## 9. Completed: `CompositeModel` removed

`CompositeModel{Nstrain,Nstress,T}` was exported but had no constructor beyond
the implicit one, no methods, and no `leafs`/`branches` fields — so every
`AbstractCompositeModel` function would have failed on it. It had zero uses
anywhere in the repository, and its docstring described it as a placeholder "for
future or experimental composite representations."

Removed:

- the name from the `export` list in [src/RheologyCalculator.jl](src/RheologyCalculator.jl)
- the struct and its docstring from [src/core/composite.jl](src/core/composite.jl)
- the Documenter entry from [docs/src/api.md](docs/src/api.md)

**Release note.** This drops a name from the public API. The package is at
v0.1.4; under Julia's 0.x convention that warrants **v0.2.0**. Bump before the
next registration, and note in the release that the removed type could not be
used with any package function.

---

## Summary

| Phase | Content | Lines removed | Risk |
|---|---|---:|---|
| 0 | Extend allocation/inference fixtures | +60 (tests) | none |
| 1 | **Fix 4 latent defects** (LTP parse, silent zero viscosity, xnorm, history bounds) | ~12 | contained |
| 2 | Delete unreachable code | ~250 | minimal |
| 3 | Fix 14 dangling refs, wrong citation, dup doc entry, shadowed arg | ~0 | none |
| 4 | Cap-family supertype, shared imports, viscosity deletions, tensor helpers public | ~215 | low |
| 5 | Deduplicate series/parallel and tensor/scalar helper pairs | ~120 | moderate |
| 6 | Replace 42 `@generated` with 3 shared tuple helpers | ~150 | moderate-high |
| 7 | `print_rheology.jl` builders, fixed buffers, `@inbounds` | ~60 | low |
| 8 | Drop 2 unused deps, enable stale-dep + docs checks, runic | — | none |
| 9 | `CompositeModel` removed | ~12 | **done** |

Net: roughly **~800 lines removed from a 4265-line `src/`** (~19%), with the
`@generated` count falling from 63 to about 11, four latent defects fixed, and
no unintended change to any computed result.

The three items with the highest value per unit of risk, if only some of this
gets done: **§1.2** (a missing viscosity method silently halves the effective
viscosity), **§4.1** (the four identical plastic model bodies), and **§2**
(250 lines that cannot run, two of which would throw if they ever did).
