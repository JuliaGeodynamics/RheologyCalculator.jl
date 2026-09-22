# Plan: assessing an implementation of the `creep.pdf` formulation

## Purpose

Decide whether the elastic history treatment in `src/` should be rebuilt on the formulation in `creep.pdf` (eqs. (1)–(6) and CASE 1/2), rather than patched in place (`ELASTIC_CORRECTION_FIX_PLAN.md`).

This plan produces **evidence and a recommendation**. It changes nothing in `src/`. Prototypes live outside `src/` and are thrown away or promoted after the decision.

References:

- `creep.pdf` (local, gitignored)
- `docs/derivations/tensor_reduction.typ`
- `docs/derivations/noncoaxial_reduction.typ`
- `docs/derivations/tensor_reduction_audit.md` (S*, D*, P*, X* IDs)
- `ELASTIC_CORRECTION_FIX_PLAN.md`

---

## 1. The candidates

| ID | Design | Summary |
|---|---|---|
| **A** | Fix in place | Keep the equation graph: block/sub-branch unknowns, spring state functions that ignore τ0, shifted local unknowns. Replace the scalar correction with one tensor correction; update history from tensors (`ELASTIC_CORRECTION_FIX_PLAN.md`). |
| **B** | PDF, linear blocks only (hybrid) | A `ParallelModel` block that holds springs and only linear elements is **collapsed**: it adds the explicit strain rate ε̇_b = (τ − Σ η*ᵢτᵒᵢ)/(2η_KV) (PDF (5) ≡ typst block equation) directly to the outer series, with **no local unknowns**. The effective strain-rate tensor ε̇* is assembled **once per step** and reduced to one invariant (PDF (6)). All other blocks keep the current graph. |
| **C** | PDF, full | Like B, plus nonlinear blocks through PDF CASE 1 (parameters depend on top-level τ_II) and CASE 2 (one extra unknown (sᵢ)_II per branch with a state-dependent element, with branch stress tensors from PDF (4)). |

What the PDF formulation implies structurally:

- Within a step, the history enters only as a **constant tensor** for linear elements. The scalar residual is r(τ_II) = τ_II − 2η* ε̇*_II.
- Branch stress tensors after the solve come from PDF (4): sᵢ = ε̇/αᵢ + (βᵢ/αᵢ)sᵢ*. These are tensors by construction, which covers the history update (S6, S7).
- There are no shifted unknowns, so viscosities are never evaluated at non-physical states (S4, S5).

---

## 2. Questions the assessment must answer

| # | Question | How it is answered |
|---|---|---|
| Q1 | Which package topologies does the PDF element cover exactly, and which need extensions? | Step 1 (mapping) |
| Q2 | Do the package's nonlinear laws fit CASE 1, CASE 2, or neither? | Step 3 |
| Q3 | What depends on the current layout of `x` (unknowns per block/sub-branch), and what breaks if blocks lose their unknowns? | Step 2 |
| Q4 | How do plastic elements and volumetric (P0) history inside blocks interact with the PDF form? | Steps 1, 3 |
| Q5 | Is B or C numerically equivalent to a full-tensor reference, including reversal and rotation? | Steps 4–6 |
| Q6 | Is B or C at least as fast as the current solver, allocation-free, type-stable, GPU-compatible, and sparsity-tracer compatible? | Steps 5–6 |
| Q7 | How much code is removed, added, or rewritten, and which public APIs break? | Step 7 |

---

## Step 1: Map the PDF element onto package composites

- [ ] Write the correspondence table, with a proof or reference for each row:

  | Package | PDF |
  |---|---|
  | `ParallelModel(SeriesModel(η₀,G₀), …, SeriesModel(η_N,G_N))` | reference + N non-reference Maxwell branches |
  | direct spring in a block | branch with η → ∞ (αᵢ = βᵢ) |
  | direct dashpot in a block | branch with G → ∞ (βᵢ = 0) |
  | several direct dashpots/springs | merged into one branch? (parallel sum; springs keep separate histories) |
  | sub-branch with >1 dashpot or >1 spring | **outside PDF** (one spring + one dashpot per branch); typst eq. for sub-branches extends it |
  | springs/dashpots in the outer series | "sequential combination" (PDF p. 4) |
  | `SeriesModel` of several blocks | sequential combination |
  | plastic element in a block (e.g. `P(DruckerPrager, η_reg)`, `test_VEVP.jl`) | **outside PDF** |
  | volumetric elastic (`Elasticity` with K, P0) in a block | **outside PDF** (deviatoric only) |
  | nesting deeper than block → sub-branch of leafs | **outside PDF** |
  | rigid element (α = 0) | degenerate in PDF (P8) |

- [ ] For each "outside PDF" row, decide: extend with the typst formulas, fall back to the graph path (hybrid B), or reject with an error.
- [ ] Record why the degenerate limits (η → ∞, G → ∞) are safe or unsafe in finite precision. Direct leafs must be dispatched explicitly, never passed as `Inf`.

Deliverable: coverage table in the assessment report (Step 8).

## Step 2: Find everything that depends on the unknown layout

Collapsing a block removes its unknowns. Every consumer of `x` needs checking.

- [ ] Inventory (from `grep x_keys|generate_equations|stress_index|branch_strain_rate_mask`):
  - `src/equation_system/`: `equations.jl`, `initial_guess.jl` (`x_keys`, `initial_guess_x`), `normalize_x.jl`, `solver.jl` (`branch_strain_rate_mask`), `tangent.jl` (`stress_index`, `tangent_block`), `inspect.jl`
  - `src/post_processing/`: `component_partition.jl`, `dissipation_partition.jl`, `post_calculations.jl`
  - tests: `test_equation_graphs.jl` (asserts graph shapes), `test_component_partition.jl`, `test_dissipation_partition.jl`, `test_initial_guess.jl`, `test_solution.jl`, `test_tangent.jl`, `test_allocations.jl`, `test_type_stability.jl`, `test_goldsby_kohlstedt.jl`, `test_nonnegative_iterate.jl`, `test_drucker_prager_dilatancy.jl`
  - docs: `composites.md`, `api.md`, `strain_rate_correction.md`
- [ ] For each consumer, classify how it relies on block unknowns: *needs them* (e.g. per-element strain-rate partitioning), *can recompute them after the solve* from PDF (4) / typst block recovery, or *indifferent*.
- [ ] Check whether component and dissipation partitions can be computed from recovered tensors rather than from `x`. That decides whether collapsing blocks is invisible to users.
- [ ] Alternative that keeps the layout: keep the block unknowns but make their equations **explicit** closures (x_b = recovered physical invariant), so `x_keys` is unchanged and those unknowns become post-processing values solved alongside. Assess its cost: Jacobian size unchanged, no speed gain.

Deliverable: an impact table (consumer → reliance → change needed → user-visible?).

## Step 3: Nonlinear semantics, plasticity, and volumetric history

- [ ] **Semantics of nonlinear laws.** Package laws are local and strain-rate-parameterized (`compute_viscosity(r; ε, …)`, `compute_strain_rate(r; τ, …)`): an element reacts to **its own** state. PDF CASE 1 makes parameters depend on the **top-level** τ_II, which is a different constitutive model (LUBBY2-style). Check whether any rheology in `src/rheology` is meant to use top-level stress. Expected result: none. If so, CASE 1 is not a drop-in; the package semantics are always CASE 2 for nonlinear elements inside blocks.
- [ ] **CASE 2 in package terms.**
  - One extra unknown per sub-branch with a state-dependent element: (sᵢ)_II ≡ t_j.
  - One per block with a state-dependent direct dashpot: e_b.
  - Law parameterization: stress (`compute_strain_rate(τ)`, closer to the PDF) or strain rate (`compute_viscosity(ε)`, implicit; typst noncoaxial §7.2 gap).
  - Choose the form that avoids a local inner solve, for each law type.
- [ ] **Plastic elements in blocks.** A block like `P(DruckerPrager, η_reg)` has a plastic multiplier unknown and a yield condition. It is not a Maxwell branch. Decide: always keep it on the graph path (hybrid), or reject history-bearing blocks that contain plasticity. Check whether any test or example combines plasticity and springs in one block (current survey: none).
- [ ] **Volumetric history.** `Elasticity` volumetric state functions already use P0 directly (`compute_volumetric_strain_rate = -(P - P0)/(K dt)`, `compute_pressure = P0 - K dt θ`). Pressure is scalar, so there is no orientation problem. Confirm that a block rewrite on the deviatoric side can leave the volumetric path untouched, and check consistency when an `Elasticity` sits in a block (deviatoric collapsed, volumetric on the graph).

Deliverable: a semantics note with a decision per law class (linear / power-law-like / plastic / volumetric).

## Step 4: Reference solutions (shared with the fix plan)

- [ ] Reuse `ELASTIC_CORRECTION_FIX_PLAN.md` Phase 0: a hand-written full-tensor time-discrete reference per topology, and load cases (constant, reversal, change of direction, 3D, scalar ε). If Phase 0 already exists, use it; otherwise build it here, once.
- [ ] Add a nonlinear reference: a monolithic Newton solve on the physical tensor unknowns (τ tensor, block strain-rate tensors, sub-branch stress tensors) with ForwardDiff, test-only. Topologies: power-law dashpot in a Maxwell sub-branch, and power-law direct dashpot in a KV block.

## Step 5: Prototype B (linear blocks) outside `src/`

Location: `examples/prototypes/pdf_blocks/` (or the gitignored benchmark folder), depending only on the package's public API plus `RheologyCalculator` internals by qualified name.

- [ ] Type-level detection of a collapsible block: it holds springs, all elements are linear (needs the `viscosity_depends_on_state` trait from fix-plan Phase 4a, prototyped locally), and sub-branches contain only leafs.
- [ ] Per step, before the Newton solve:
  1. coefficients η_M,j, η*ᵢ, η_KV,b (all constant for linear elements)
  2. ε̇* tensor = ε̇ + H, reduced once to ε̇*_II
  3. the collapsed block contributes its viscous part τ/(2η_KV,b) to the outer series residual
- [ ] Scalar residual in τ_II (plus remaining graph unknowns for non-collapsed blocks). Solve it with the existing Newton / line search if possible.
- [ ] After the solve: τ = τ_II n, block tensors, spring tensors (PDF (4) / typst §6).
- [ ] Compare with the Step 4 reference at every step (τ and every spring tensor), and with the current `solve` for coaxial constant loading (should match to round-off).
- [ ] Measure with `BenchmarkTools`: `solve` time and iterations for Maxwell, KV, mixed, Burgers and two-block models, against the current code; allocations (`@allocated` == 0 after warm-up); `@inferred`; a ForwardDiff Jacobian of the reduced residual; `SparseConnectivityTracer` on the reduced residual.
- [ ] Record the Jacobian size before and after (unknowns removed per model).

## Step 6: Prototype the CASE 2 subset (only if Step 3 says it is needed)

- [ ] One topology: `S(η₁, P(η₂, S(PowerLaw, G)))`. Unknowns τ_II and t_j. Residuals as in PDF CASE 2 step 4 / typst `rx0`, `rxj`, with tensors assembled at each residual evaluation.
- [ ] Compare with the nonlinear reference (Step 4) under constant loading, reversal, and a change of direction.
- [ ] The same measurements as Step 5, plus Newton iteration counts against the current graph solver on the same model (which is known to be wrong for S4; compare cost, not answers).

## Step 7: Integration cost estimate

For B, and for C if Step 6 ran:

- [ ] List the `src/` functions that would be removed, rewritten, or added, with rough line counts. Candidates for removal: `_kv_implicit_*`, `subtract_elastic_correction`, `_direct_leaf_elastic_correction`, `_branch_elastic_stress`, `_branch_elastic_info`, and the shifted-unknown handling in `compute_stress_elastic`.
- [ ] List public API changes: `x_keys` / `initial_guess_x` / `normalisation_x` for models with collapsed blocks, `effective_strain_rate_correction`, `elastic_stress_history_*`, `compute_stress_elastic`, `tangent*`, `inspect`, the partitions.
- [ ] Estimate the test churn (`test_equation_graphs.jl` expectations etc.) and the docs rewrite (`composites.md`, `strain_rate_correction.md`).
- [ ] Identify what B and C make **unnecessary** from the fix plan (for example, with B the correction is constant per step for linear blocks, so the x-dependent H(x) of fix-plan Phase 2 becomes obsolete).

## Step 8: Report and decision

Write `docs/derivations/pdf_implementation_assessment.md` (local until reviewed) with:

- [ ] the coverage table (Step 1), impact table (Step 2), semantics decisions (Step 3)
- [ ] the numerical equivalence results (Steps 5–6): max relative error per topology and load case
- [ ] performance and constraint results (Steps 5–6): time, iterations, allocations, inference, GPU/tracer compatibility
- [ ] the integration cost (Step 7)
- [ ] a decision matrix:

  | Criterion | A (fix in place) | B (PDF, linear blocks) | C (PDF, full) |
  |---|---|---|---|
  | Correct for linear, non-coaxial, reversal | | | |
  | Correct for nonlinear in blocks | | | |
  | Topology coverage (incl. plastic, volumetric) | | | |
  | Solve cost (time, Jacobian size) | | | |
  | Allocation/inference/GPU/tracer | | | |
  | Public API breakage | | | |
  | `src/` lines removed / added | | | |
  | Risk / effort | | | |

- [ ] A recommendation, and the ordered next steps (which may merge parts of `ELASTIC_CORRECTION_FIX_PLAN.md`).

---

## Stop conditions

End the assessment early and report if any of these holds:

- Step 2 shows that removing block unknowns breaks a user-facing feature (for example component or dissipation partitioning) that cannot be recovered after the solve. Then B/C only make sense with the layout-preserving variant, and the speed argument disappears.
- Step 5 cannot match the reference for linear reversal or rotation to round-off. That would mean the formulation mapping has an error.
- Step 3 finds that package semantics need CASE 2 for common models and that the extra unknowns are equivalent to today's block unknowns. Then C converges back to A with physical instead of shifted unknowns, which is a variant of A, not a new design.

## Decisions needed before starting

| ID | Question | Default |
|---|---|---|
| R1 | Where do prototypes live (`examples/prototypes/`, gitignored benchmark folder, or a branch)? | a local branch, `examples/prototypes/pdf_blocks/` |
| R2 | Is changing the unknown layout for models with springs in blocks acceptable in principle (breaking)? | assess, don't assume |
| R3 | Scope of nonlinear support to assess: none (B only), or CASE 2 for power-law-like laws (C)? | B first; C only if Step 3 finds a real use case |
